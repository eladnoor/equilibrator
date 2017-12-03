import logging
import numpy
import pulp

from matching import query_parser
from util import constants

DEFAULT_RT = constants.R * constants.DEFAULT_TEMP


class MDFResult(object):

    def __init__(self, model, mdf,
                 concentrations, dG0_r_cov_eigen,
                 reaction_prices, compound_prices,
                 min_total_dG=None, max_total_dG=None):
        """Initialize.

        Args:
            model: PathwayModel generating these results.
            mdf: value of the result.
            concentrations: metabolite concentrations at optimum.
            dG0_r_cov_eigen: dGr0' covariance eigenvalues at optimum.
            reaction_prices: shadow prices for reactions.
            compound_prices: shadow prices for compounds.
        """
        self.model = model
        self.mdf = mdf
        self.concentrations = concentrations
        self.dG0_r_cov_eigen = dG0_r_cov_eigen
        self.reaction_prices = reaction_prices
        self.compound_prices = compound_prices

        self.dG_r_prime = model.CalculateReactionEnergiesUsingConcentrations(
            concentrations)
        self.dG_r_prime_raw = self.dG_r_prime + numpy.dot(model.dG0_r_std,
                                                       dG0_r_cov_eigen)

        # adjust dG to flux directions
        self.dG_r_prime_adj = model.I_dir * self.dG_r_prime_raw

        # May be set after initialization. Optional.
        self.max_total_dG = max_total_dG
        self.min_total_dG = min_total_dG


class PathwayThermoModel(object):
    """Container for doing pathway-level thermodynamic analysis."""

    DEFAULT_FORMATION_LB = -1e6
    DEFAULT_FORMATION_UB = 1e6
    DEFAULT_REACTION_LB = -1e3
    DEFAULT_REACTION_UB = 0.0
    DEFAULT_C_RANGE = (1e-6, 0.1)
    DEFAULT_PHYSIOLOGICAL_CONC = 1e-3

    def __init__(self, S, fluxes, dG0_r_prime, cids, rids,
                 dG0_r_std=None, concentration_bounds=None):
        """Create a pathway object.

        Args:
            S: Stoichiometric matrix of the pathway.
                Reactions are on the rows, compounds on the columns.
            fluxes: the list of relative fluxes through each of the reactions.
                By default, all fluxes are 1.
            dG0_r_prime: the change in Gibbs energy for the reactions
                in standard conditions, corrected for pH, ionic strength, etc.
                Should be a column vector in numpy.matrix format.
            dG0_r_std: (optional) the square root of the covariance matrix
                corresponding to the uncertainty in the dG0_r values.
            concentration_bounds: a bounds.Bounds object expressing bounds on
                the metabolite concentrations.
        """
        self.pulp_solver = pulp.GLPK_CMD(msg=0, options=['--xcheck'])

        self.S = S
        self.Nc, self.Nr = S.shape

        self.dG0_r_prime = dG0_r_prime
        if dG0_r_std is None:
            self.dG0_r_std = numpy.matrix(numpy.zeros((self.Nr, self.Nr)))
        else:
            self.dG0_r_std = dG0_r_std

        # Make sure dG0_r' is the right size
        assert self.dG0_r_prime.shape[0] == self.Nr
        assert self.dG0_r_std.shape[0] == self.Nr
        assert self.dG0_r_std.shape[1] == self.Nr

        self.fluxes = fluxes
        if self.fluxes is None:
            self.fluxes = numpy.ones((1, self.Nr))
        self.fluxes = numpy.matrix(self.fluxes)
        assert self.fluxes.shape[1] == self.Nr, 'Fluxes required for all reactions'

        self.I_dir = numpy.matrix(numpy.diag(map(numpy.sign, self.fluxes.flat)))
        self.Nr_active = int(sum(self.fluxes.T != 0))

        self.cids = cids
        self.rids = rids
        self.concentration_bounds = concentration_bounds

        if self.concentration_bounds is None:
            lb, ub = self.DEFAULT_C_RANGE
            self.concentration_bounds = bounds.Bounds(default_lb=lb, default_ub=ub)

        # Currently unused bounds on reaction dGs.
        self.r_bounds = None

    def CalculateReactionEnergiesUsingConcentrations(self, concentrations):
        log_conc = numpy.log(concentrations)
        if numpy.isnan(self.dG0_r_prime).any():
            dG_r_prime = self.dG0_r_prime.copy()
            for r in range(self.Nr):
                reactants = list(self.S[:, r].nonzero()[0].flat)
                dG_r_prime[0, r] += DEFAULT_RT * log_conc[reactants, 0].T * self.S[reactants, r]
            return dG_r_prime
        else:
            return self.dG0_r_prime + DEFAULT_RT * self.S.T * log_conc

    def GetPhysiologicalConcentrations(self, bounds=None):
        conc = numpy.matrix(numpy.ones((self.Nc, 1))) * self.DEFAULT_PHYSIOLOGICAL_CONC
        if bounds:
            for i, bound in enumerate(bounds):
                lb, ub = bound
                if lb is not None and ub is not None:
                    if not (lb < conc[i, 0] < ub):
                        conc[i, 0] = numpy.sqrt(lb * ub)

        return conc

    def _MakeLnConcentratonBounds(self):
        """Make bounds on logarithmic concentrations.

        Returns:
            A two-tuple (lower bounds, upper bounds).
        """
        bounds = self.concentration_bounds.GetLnBounds(self.cids)
        return bounds

    def _MakeDrivingForceConstraints(self, ln_conc_lb, ln_conc_ub):
        """Generates the A matrix and b & c vectors that can be used in a
        standard form linear problem:
                max          c'x
                subject to   Ax <= b

        x is the vector of (y | log-conc | B)
        where y dG'0 are the reaction Gibbs energy variables, log-conc
        are the natural log of the concentrations of metabolites, and
        B is the max-min driving force variable which is being maximized
        by the LP
        """
        inds = numpy.nonzero(numpy.diag(self.I_dir))[0].tolist()

        # driving force
        A11 = self.I_dir[inds] * self.dG0_r_std
        A12 = self.I_dir[inds] * self.S.T * DEFAULT_RT
        A13 = numpy.ones((len(inds), 1))

        # covariance var ub and lb
        A21 = numpy.eye(self.Nr)
        A22 = numpy.zeros((self.Nr, self.Nc))
        A23 = numpy.zeros((self.Nr, 1))

        # log conc ub and lb
        A31 = numpy.zeros((self.Nc, self.Nr))
        A32 = numpy.eye(self.Nc)
        A33 = numpy.zeros((self.Nc, 1))

        # upper bound values
        b1 = -self.I_dir[inds] * self.dG0_r_prime
        b2 = numpy.ones((self.Nr, 1))

        A = numpy.matrix(numpy.vstack([numpy.hstack([ A11,  A12,  A13]),   # driving force
                                 numpy.hstack([ A21,  A22,  A23]),   # covariance var ub
                                 numpy.hstack([-A21,  A22,  A23]),   # covariance var lb
                                 numpy.hstack([ A31,  A32,  A33]),   # log conc ub
                                 numpy.hstack([ A31, -A32,  A33])])) # log conc lb

        b = numpy.matrix(numpy.vstack([b1, b2, b2, ln_conc_ub, -ln_conc_lb]))

        c = numpy.matrix(numpy.zeros((A.shape[1], 1)))
        c[-1, 0] = 1.0

        # change the constaints such that reaction that have an explicit
        # r_bound will not be constrained by B, but will be constained by
        # their specific bounds. Note that we need to divide the bound
        # by R*T since the variables in the LP are not in kJ/mol but in units
        # of R*T.
        if self.r_bounds:
            for i, r_ub in enumerate(self.r_bounds):
                if r_ub is not None:
                    A[i, -1] = 0.0
                    b[i, 0] += r_ub

        return A, b, c

    def _GetPrimalVariablesAndConstants(self):
        # Define and apply the constraints on the concentrations
        ln_conc_lb, ln_conc_ub = self._MakeLnConcentratonBounds()

        # Create the driving force variable and add the relevant constraints
        A, b, c = self._MakeDrivingForceConstraints(ln_conc_lb, ln_conc_ub)

        # the dG'0 covariance eigenvariables
        y = pulp.LpVariable.dicts("y", ["%d" % i for i in range(self.Nr)])
        y = [y["%d" % i] for i in range(self.Nr)]

        # ln-concentration variables
        l = pulp.LpVariable.dicts("l", ["%d" % i for i in range(self.Nc)])
        l = [l["%d" % i] for i in range(self.Nc)]

        return A, b, c, y, l

    def _GetDualVariablesAndConstants(self):
        # Define and apply the constraints on the concentrations
        ln_conc_lb, ln_conc_ub = self._MakeLnConcentratonBounds()

        # Create the driving force variable and add the relevant constraints
        A, b, c = self._MakeDrivingForceConstraints(ln_conc_lb, ln_conc_ub)

        w = pulp.LpVariable.dicts("w",
                                  ["%d" % i for i in range(self.Nr_active)],
                                  lowBound=0)
        w = [w["%d" % i] for i in range(self.Nr_active)]

        g = pulp.LpVariable.dicts("g",
                                  ["%d" % i for i in range(2*self.Nr)],
                                  lowBound=0)
        g = [g["%d" % i] for i in range(2*self.Nr)]

        z = pulp.LpVariable.dicts("z",
                                  ["%d" % i for i in range(self.Nc)],
                                  lowBound=0)
        z = [z["%d" % i] for i in range(self.Nc)]

        u = pulp.LpVariable.dicts("u",
                                  ["%d" % i for i in range(self.Nc)],
                                  lowBound=0)
        u = [u["%d" % i] for i in range(self.Nc)]

        return A, b, c, w, g, z, u

    def _GetTotalEnergyProblem(self,
                               min_driving_force=0.0,
                               objective=pulp.LpMinimize):

        A, b, _c, y, l = self._GetPrimalVariablesAndConstants()
        x = y + l + [min_driving_force]
        lp = pulp.LpProblem("MDF", objective)

        for j in range(3*self.Nr + 2 * self.Nc):
            row = [A[j, i] * x[i] for i in range(self.Nr + self.Nc + 1)]
            lp += (pulp.lpSum(row) <= b[j, 0]), "energy_%02d" % j

        total_g = pulp.LpVariable("g_tot")
        total_g0 = float(self.fluxes * self.dG0_r_prime)
        total_reaction = self.S * self.fluxes.T
        row = [total_reaction[i, 0] * x[i] for i in range(self.Nc)]
        lp += (total_g == total_g0 + pulp.lpSum(row)), "Total G"

        lp.setObjective(total_g)

        return lp, total_g

    def _MakeMDFProblem(self):
        """Create a CVXOPT problem for finding the Maximal Thermodynamic
        Driving Force (MDF).

        Does not set the objective function... leaves that to the caller.

        Returns:
            the linear problem object, and the three types of variables as arrays
        """
        A, b, c, y, l = self._GetPrimalVariablesAndConstants()
        B = pulp.LpVariable("mdf")
        x = y + l + [B]
        lp = pulp.LpProblem("MDF_PRIMAL", pulp.LpMaximize)

        cnstr_names = ["driving_force_%02d" % j for j in range(self.Nr_active)] + \
                      ["covariance_var_ub_%02d" % j for j in range(self.Nr)] + \
                      ["covariance_var_lb_%02d" % j for j in range(self.Nr)] + \
                      ["log_conc_ub_%02d" % j for j in range(self.Nc)] + \
                      ["log_conc_lb_%02d" % j for j in range(self.Nc)]

        for j in range(A.shape[0]):
            row = [A[j, i] * x[i] for i in range(A.shape[1])]
            lp += (pulp.lpSum(row) <= b[j, 0]), cnstr_names[j]

        objective = pulp.lpSum([c[i] * x[i] for i in range(A.shape[1])])
        lp.setObjective(objective)

        return lp, objective, y, l, B

    def _MakeMDFProblemDual(self):
        """Create a CVXOPT problem for finding the Maximal Thermodynamic
        Driving Force (MDF).

        Does not set the objective function... leaves that to the caller.

        Returns:
            the linear problem object, and the four types of variables as arrays
        """
        A, b, c, w, g, z, u = self._GetDualVariablesAndConstants()
        x = w + g + z + u
        lp = pulp.LpProblem("MDF_DUAL", pulp.LpMinimize)

        cnstr_names = ["y_%02d" % j for j in range(self.Nr)] + \
                      ["l_%02d" % j for j in range(self.Nc)] + \
                      ["MDF"]

        for i in range(A.shape[1]):
            row = [A[j, i] * x[j] for j in range(A.shape[0])]
            lp += (pulp.lpSum(row) == c[i, 0]), cnstr_names[i]

        objective = pulp.lpSum([b[i] * x[i] for i in range(A.shape[0])])
        lp.setObjective(objective)

        return lp, objective, w, g, z, u

    def FindMDF(self, calculate_totals=True):
        """Find the MDF (Optimized Bottleneck Energetics).

        Args:
            c_range: a tuple (min, max) for concentrations (in M).
            bounds: a list of (lower bound, upper bound) tuples for compound
                concentrations.

        Returns:
            A 2 (optimal dGfs, optimal concentrations, optimal mdf).
        """
        lp_primal, primal_obj, y, l, B = self._MakeMDFProblem()
        lp_primal.solve(self.pulp_solver)

        if lp_primal.status != pulp.LpStatusOptimal:
            logging.warning('LP status %s', lp_primal.status)
            raise pulp.solvers.PulpSolverError(
                "Cannot solve MDF primal optimization problem")

        y = numpy.matrix(map(pulp.value, y)).T
        l = numpy.matrix(map(pulp.value, l)).T
        mdf = pulp.value(B)
        conc = numpy.exp(l)
        dG0_r_prime = self.dG0_r_prime + numpy.dot(self.dG0_r_std, y)

        lp_dual, dual_obj, w, g, z, u = self._MakeMDFProblemDual()
        lp_dual.solve(self.pulp_solver)
        if lp_dual.status != pulp.LpStatusOptimal:
            raise pulp.solvers.PulpSolverError("cannot solve MDF dual")

        primal_obj = pulp.value(primal_obj)
        dual_obj = pulp.value(dual_obj)
        if abs(primal_obj - dual_obj) > 1e-3:
            raise pulp.solvers.PulpSolverError("Primal != Dual (%.5f != %.5f)"
            % (primal_obj, dual_obj))

        w = map(pulp.value, w)
        z = map(pulp.value, z)
        u = map(pulp.value, u)
        reaction_prices = numpy.matrix(w).T
        compound_prices = numpy.matrix(z).T - numpy.matrix(u).T

        ret = MDFResult(self, mdf, conc, y, reaction_prices, compound_prices)

        if calculate_totals:
            # find the maximum and minimum total Gibbs energy of the pathway,
            # under the constraint that the driving force of each reaction is >= MDF
            lp_total, total_dg = self._GetTotalEnergyProblem(mdf - 1e-2, pulp.LpMinimize)
            lp_total.solve(self.pulp_solver)
            if lp_total.status != pulp.LpStatusOptimal:
                logging.warning("cannot solve minimal total delta-G problem")
                ret.min_total_dG = numpy.nan
            else:
                ret.min_total_dG = pulp.value(total_dg)

            lp_total, total_dg = self._GetTotalEnergyProblem(mdf - 1e-2, pulp.LpMaximize)
            lp_total.solve(self.pulp_solver)
            if lp_total.status != pulp.LpStatusOptimal:
                logging.warning("cannot solve maximal total delta-G problem")
                ret.max_total_dG = numpy.nan
            else:
                ret.max_total_dG = pulp.value(total_dg)

        return ret

    @property
    def mdf_result(self):
        ret = self.FindMDF()
        return ret
