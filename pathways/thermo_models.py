import csv
import logging
import numpy as np
import pulp
import re

from matching import query_parser
from gibbs import constants
from gibbs.reaction import Reaction, CompoundWithCoeff

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
        self.dG_r_prime_raw = self.dG_r_prime + np.dot(model.dG0_r_std,
                                                       dG0_r_cov_eigen)

        # adjust dG to flux directions
        self.dG_r_prime_adj = model.I_dir * self.dG_r_prime_raw
        
        # May be set after initialization. Optional.
        self.max_total_dG = max_total_dG
        self.min_total_dG = min_total_dG

    @property
    def reaction_dGs_and_prices(self):
        return zip(self.dG_r_prime_adj.flatten().tolist()[0],
                   self.reaction_prices.flatten().tolist()[0])
    

class PathwayModel(object):
    """Container for doing pathway-level thermodynamic analysis."""
   
    DEFAULT_FORMATION_LB = -1e6
    DEFAULT_FORMATION_UB = 1e6
    DEFAULT_REACTION_LB = -1e3
    DEFAULT_REACTION_UB = 0.0
    DEFAULT_C_RANGE = (1e-6, 0.1)
    DEFAULT_PHYSIOLOGICAL_CONC = 1e-3
   
    def __init__(self, S, dG0_r_prime, dG0_r_std=None, fluxes=None):
        """Create a pathway object.
       
        Args:
            S: Stoichiometric matrix of the pathway.
                Reactions are on the rows, compounds on the columns.
            dG0_r_prime: the change in Gibbs energy for the reactions
                in standard conditions, corrected for pH, ionic strength, etc.
                Should be a column vector in numpy.matrix format.
            dG0_r_std: (optional) the square root of the covariance matrix
                corresponding to the uncertainty in the dG0_r values.
            fluxes: the list of relative fluxes through each of the reactions.
                By default, all fluxes are 1.
        """
        self.pulp_solver = pulp.GLPK_CMD(msg=0, options=['--xcheck'])
        
        self.S = S
        self.Nc, self.Nr = S.shape

        self.dG0_r_prime = dG0_r_prime
        if dG0_r_std is None:
            self.dG0_r_std = np.matrix(np.zeros((self.Nr, self.Nr)))
        else:
            self.dG0_r_std = dG0_r_std
       
        # Make sure dG0_r' is the right size
        assert self.dG0_r_prime.shape[0] == self.Nr
        assert self.dG0_r_std.shape[0] == self.Nr
        assert self.dG0_r_std.shape[1] == self.Nr

        self.fluxes = fluxes
        if self.fluxes is None:
            self.fluxes = np.ones((1, self.Nr))
        self.fluxes = np.matrix(self.fluxes)
        assert self.fluxes.shape[1] == self.Nr, 'Fluxes required for all reactions'
        
        self.I_dir = np.matrix(np.diag(map(np.sign, self.fluxes.flat)))
        self.Nr_active = int(sum(self.fluxes.T != 0))
        self.c_bounds = None
        self.r_bounds = None
        self.c_range = self.DEFAULT_C_RANGE

    def CalculateReactionEnergiesUsingConcentrations(self, concentrations):
        log_conc = np.log(concentrations)
        if np.isnan(self.dG0_r_prime).any():
            dG_r_prime = self.dG0_r_prime.copy()
            for r in xrange(self.Nr):
                reactants = list(self.S[:, r].nonzero()[0].flat)
                dG_r_prime[0, r] += DEFAULT_RT * log_conc[reactants, 0].T * self.S[reactants, r]
            return dG_r_prime
        else:
            return self.dG0_r_prime + DEFAULT_RT * self.S.T * log_conc

    def GetPhysiologicalConcentrations(self, bounds=None):
        conc = np.matrix(np.ones((self.Nc, 1))) * self.DEFAULT_PHYSIOLOGICAL_CONC
        if bounds:
            for i, bound in enumerate(bounds):
                lb, ub = bound
                if lb is not None and ub is not None:
                    if not (lb < conc[i, 0] < ub):
                        conc[i, 0] = np.sqrt(lb * ub)
       
        return conc

    def _MakeLnConcentratonBounds(self):
        """Make bounds on logarithmic concentrations.
        
        Returns:
            A two-tuple (lower bounds, upper bounds).
        """
        c_lower, c_upper = self.c_range or self.DEFAULT_C_RANGE
        ln_conc_lb = np.matrix(np.ones((self.Nc, 1)) * np.log(c_lower))
        ln_conc_ub = np.matrix(np.ones((self.Nc, 1)) * np.log(c_upper))
       
        if self.c_bounds:
            for i, bound in enumerate(self.c_bounds):
                lb, ub = bound
                log_lb = np.log(lb or c_lower)
                log_ub = np.log(ub or c_upper)
                if log_lb > log_ub:
                    raise Exception("Lower bound is greater than upper bound: "
                                    "%d > %d" % (log_lb, log_ub))
                elif abs(log_lb - log_ub) < 1e-2:
                    log_lb = log_ub - 1e-2
                   
                ln_conc_lb[i, 0] = log_lb
                ln_conc_ub[i, 0] = log_ub

        return ln_conc_lb, ln_conc_ub

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
        inds = np.nonzero(np.diag(self.I_dir))[0].tolist()
        
        # driving force
        A11 = self.I_dir[inds] * self.dG0_r_std
        A12 = self.I_dir[inds] * self.S.T * DEFAULT_RT
        A13 = np.ones((len(inds), 1))
        
        # covariance var ub and lb
        A21 = np.eye(self.Nr)
        A22 = np.zeros((self.Nr, self.Nc))
        A23 = np.zeros((self.Nr, 1))
        
        # log conc ub and lb
        A31 = np.zeros((self.Nc, self.Nr))
        A32 = np.eye(self.Nc)
        A33 = np.zeros((self.Nc, 1))
        
        # upper bound values
        b1 = -self.I_dir[inds] * self.dG0_r_prime
        b2 = np.ones((self.Nr, 1))
        
        A = np.matrix(np.vstack([np.hstack([ A11,  A12,  A13]),   # driving force
                                 np.hstack([ A21,  A22,  A23]),   # covariance var ub 
                                 np.hstack([-A21,  A22,  A23]),   # covariance var lb 
                                 np.hstack([ A31,  A32,  A33]),   # log conc ub
                                 np.hstack([ A31, -A32,  A33])])) # log conc lb

        b = np.matrix(np.vstack([b1, b2, b2, ln_conc_ub, -ln_conc_lb]))

        c = np.matrix(np.zeros((A.shape[1], 1)))
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
        y = pulp.LpVariable.dicts("y", ["%d" % i for i in xrange(self.Nr)])
        y = [y["%d" % i] for i in xrange(self.Nr)]

        # ln-concentration variables
        l = pulp.LpVariable.dicts("l", ["%d" % i for i in xrange(self.Nc)])
        l = [l["%d" % i] for i in xrange(self.Nc)]

        return A, b, c, y, l
   
    def _GetDualVariablesAndConstants(self):
        # Define and apply the constraints on the concentrations
        ln_conc_lb, ln_conc_ub = self._MakeLnConcentratonBounds()

        # Create the driving force variable and add the relevant constraints
        A, b, c = self._MakeDrivingForceConstraints(ln_conc_lb, ln_conc_ub)
       
        w = pulp.LpVariable.dicts("w", 
                                  ["%d" % i for i in xrange(self.Nr_active)],
                                  lowBound=0)
        w = [w["%d" % i] for i in xrange(self.Nr_active)]

        g = pulp.LpVariable.dicts("g", 
                                  ["%d" % i for i in xrange(2*self.Nr)],
                                  lowBound=0)
        g = [g["%d" % i] for i in xrange(2*self.Nr)]

        z = pulp.LpVariable.dicts("z", 
                                  ["%d" % i for i in xrange(self.Nc)],
                                  lowBound=0)
        z = [z["%d" % i] for i in xrange(self.Nc)]

        u = pulp.LpVariable.dicts("u", 
                                  ["%d" % i for i in xrange(self.Nc)],
                                  lowBound=0)
        u = [u["%d" % i] for i in xrange(self.Nc)]
        
        return A, b, c, w, g, z, u
   
    def _GetTotalEnergyProblem(self,
                               min_driving_force=0.0,
                               objective=pulp.LpMinimize):
        
        A, b, _c, y, l = self._GetPrimalVariablesAndConstants()
        x = y + l + [min_driving_force]
        lp = pulp.LpProblem("MDF", objective)
        
        for j in xrange(3*self.Nr + 2 * self.Nc):
            row = [A[j, i] * x[i] for i in xrange(self.Nr + self.Nc + 1)]
            lp += (pulp.lpSum(row) <= b[j, 0]), "energy_%02d" % j
        
        total_g = pulp.LpVariable("g_tot")
        total_g0 = float(self.fluxes * self.dG0_r_prime)
        total_reaction = self.S * self.fluxes.T
        row = [total_reaction[i, 0] * x[i] for i in xrange(self.Nc)]
        lp += (total_g == total_g0 + pulp.lpSum(row)), "Total G"

        lp.setObjective(total_g)
        
        #lp.writeLP("res/total_g.lp")
        
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
        
        cnstr_names = ["driving_force_%02d" % j for j in xrange(self.Nr_active)] + \
                      ["covariance_var_ub_%02d" % j for j in xrange(self.Nr)] + \
                      ["covariance_var_lb_%02d" % j for j in xrange(self.Nr)] + \
                      ["log_conc_ub_%02d" % j for j in xrange(self.Nc)] + \
                      ["log_conc_lb_%02d" % j for j in xrange(self.Nc)]
          
        for j in xrange(A.shape[0]):
            row = [A[j, i] * x[i] for i in xrange(A.shape[1])]
            lp += (pulp.lpSum(row) <= b[j, 0]), cnstr_names[j]
        
        objective = pulp.lpSum([c[i] * x[i] for i in xrange(A.shape[1])])
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

        cnstr_names = ["y_%02d" % j for j in xrange(self.Nr)] + \
                      ["l_%02d" % j for j in xrange(self.Nc)] + \
                      ["MDF"]
        
        for i in xrange(A.shape[1]):
            row = [A[j, i] * x[j] for j in xrange(A.shape[0])]
            lp += (pulp.lpSum(row) == c[i, 0]), cnstr_names[i]

        objective = pulp.lpSum([b[i] * x[i] for i in xrange(A.shape[0])])
        lp.setObjective(objective)
        
        #lp.writeLP("res/mdf_dual.lp")
        
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
            raise pulp.solvers.PulpSolverError("cannot solve MDF primal")
        y = np.matrix(map(pulp.value, y)).T
        l = np.matrix(map(pulp.value, l)).T
        mdf = pulp.value(B)
        conc = np.exp(l)
        dG0_r_prime = self.dG0_r_prime + np.dot(self.dG0_r_std, y)

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
        reaction_prices = np.matrix(w).T
        compound_prices = np.matrix(z).T - np.matrix(u).T
                  
        ret = MDFResult(self, mdf, conc, y, reaction_prices, compound_prices)

        if calculate_totals:
            # find the maximum and minimum total Gibbs energy of the pathway,
            # under the constraint that the driving force of each reaction is >= MDF
            lp_total, total_dg = self._GetTotalEnergyProblem(mdf - 1e-2, pulp.LpMinimize)
            lp_total.solve(self.pulp_solver)
            if lp_total.status != pulp.LpStatusOptimal:
                logging.warning("cannot solve minimal total delta-G problem")
                ret.min_total_dG = np.nan
            else:
                ret.min_total_dG = pulp.value(total_dg)
        
            lp_total, total_dg = self._GetTotalEnergyProblem(mdf - 1e-2, pulp.LpMaximize)
            lp_total.solve(self.pulp_solver)
            if lp_total.status != pulp.LpStatusOptimal:
                logging.warning("cannot solve maximal total delta-G problem")
                ret.max_total_dG = np.nan
            else:
                ret.max_total_dG = pulp.value(total_dg)
        
        return ret
    
    @property
    def mdf_result(self):
        ret = self.FindMDF()
        return ret


class KeggPathwayModel(PathwayModel):
    """TODO: why is this separate from PathwayModel? Also, constructor
    arguments are in a strange order.
    """

    def __init__(self, S, dG0_r_prime, cids, rids,
                 dG0_r_std=None,
                 fluxes=None,
                 rid2bounds=None,
                 cid2bounds=None,
                 c_range=None):
        """
            S           - the stoichiometric matrix
            rid         - a list of names of the reactions in S
            fluxes      - a vector the relative fluxes in each reaction
            cids        - a list of names of the compounds in S
            formation_energies - the standard Gibbs energy of formation of the 
                                 compounds
            rid2bounds - a dictionary mapping rid to an upper bound on its dG'
                         if the value is None then the upper bound
                         is the B variable (corresponding to the MDF)
            reaction_energies - the standard Gibbs energies of the reactions
            cid2bounds - a dictionary mapping cid to a pair of lower/upper
                         bound on its concentration. if the value is (None, None)
                         the default bounds are used (i.e. c_range)
            c_range    - the default lower/upper bounds on the compound conc.
        """
        PathwayModel.__init__(self, S, dG0_r_prime=dG0_r_prime, dG0_r_std=dG0_r_std, fluxes=fluxes)
        assert len(cids) == self.Nc
        assert len(rids) == self.Nr
       
        self.rids = rids
        self.cids = cids
        if cid2bounds:
            self.c_bounds = [cid2bounds.get(cid, (None, None)) for cid in self.cids]
        else:
            self.c_bounds = None
        self.cid2bounds = cid2bounds or {}
        
        if rid2bounds:
            self.r_bounds = [rid2bounds.get(rid, None) for rid in self.rids]
        else:
            self.r_bounds = None
        
        self.rid2bounds = rid2bounds or {}
        self.c_range = c_range or self.DEFAULT_C_RANGE

    def GetConcentrationBounds(self, cid):
        lb, ub = self.c_range
        if cid in self.cid2bounds:
            lb, ub = self.cid2bounds[cid]
        if lb is None:
            lb = self.c_range[0]
        if ub is None:
            ub = self.c_range[1]
        return lb, ub