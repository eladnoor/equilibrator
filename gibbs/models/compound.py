#!/usr/bin/python

import json
import logging
import numpy
import re
import base64

from django.http import Http404
from django.db import models
from django.utils.text import slugify
from django.apps import apps
from util import constants
from .. import formula_parser
from .. import conditions
from util.thumbnail import InChI2Thumbnail


class CommonName(models.Model):
    """
        A common name of a compound.
    """
    name = models.CharField(max_length=500)
    enabled = models.BooleanField(default=True)

    @staticmethod
    def GetOrCreate(name):
        """Attempts to fetch the CommonName object, or creates it if not present.

        Args:
            name: the name to use.

        Returns:
            A CommonName object.
        """
        try:
            n = CommonName.objects.get(name=name, enabled=True)
            return n
        except Exception:
            n = CommonName(name=name, enabled=True)
            n.save()
            return n

    def __unicode__(self):
        return self.name

    def TypeStr(self):
        if self.compound_set.count():
            return 'Compound'
        elif self.enzyme_set.count():
            return 'Enzyme'
        return ''


class ValueSource(models.Model):
    """
        The source of a particular numeric value.
    """
    # A JSON of all the relevant data.
    name = models.CharField(max_length=30)

    # A JSON of all the relevant data.
    data = models.CharField(max_length=2048)

    def __init__(self, *args, **kwargs):
        super(ValueSource, self).__init__(*args, **kwargs)
        self._data_dict = json.loads(self.data)

    def __unicode__(self):
        return self.name

    def __str__(self):
        return self.name

    def __eq__(self, other):
        """Equality operator."""
        if not other:
            return False

        if not hasattr(other, 'name'):
            return False

        return self.name == other.name

    def GetAuthorList(self):
        return ', '.join(self._data_dict.get('author', []))
    author = property(GetAuthorList)

    def GetYear(self):
        return self._data_dict.get('year', '2000 BC')
    year = property(GetYear)

    def GetTitle(self):
        return self._data_dict.get('title', None)
    title = property(GetTitle)

    def GetJournal(self):
        return self._data_dict.get('journal', '')
    journal = property(GetJournal)

    def GetURL(self):
        return self._data_dict.get('url', None)
    url = property(GetURL)

    def GetDOI(self):
        return self._data_dict.get('doi', None)
    doi = property(GetDOI)

    def GetPMID(self):
        return self._data_dict.get('pmid', None)
    pmid = property(GetPMID)

    def GetCitation(self):
        if not self.title:
            return None
        if self.url is not None:
            return '%s. <a href="%s"><strong>%s</strong></a>. <i>%s</i>' % \
                (self.author, self.url, self.title, self.journal)
        else:
            return '%s. <strong>%s</strong>. <i>%s</i>' % \
                (self.author, self.title, self.journal)
    citation = property(GetCitation)


class Specie(models.Model):
    """
        A single specie of a compound.
    """
    # The ID of the compound in KEGG.
    kegg_id = models.CharField(max_length=10)

    # The number of hydrogens in the species.
    number_of_hydrogens = models.IntegerField()

    # The number of Mg2+ ions bound to this species.
    number_of_mgs = models.IntegerField(default=0)

    # The net charge (eV).
    net_charge = models.IntegerField()

    # The standard formation energy in kJ/mol.
    formation_energy = models.FloatField()

    # The source of this value.
    formation_energy_source = models.ForeignKey(ValueSource,
                                                on_delete=models.CASCADE)

    # The phase (s, l, g, or aq)
    phase = models.CharField(max_length=10,
                             choices=constants.PHASE_CHOICES,
                             default=constants.DEFAULT_PHASE)

    def Transform(self, aq_params):
        """
            Transform this individual estimate to difference conditions.
        """

        sqrt_I = numpy.sqrt(aq_params.ionic_strength)
        pH = aq_params.pH
        pMg = aq_params.pMg
        nH = self.number_of_hydrogens
        nMg = self.number_of_mgs
        z = self.net_charge
        dG = self.formation_energy
        dG_prime = dG

        # add the potential related to the pH
        if nH > 0:
            dG_prime += nH * constants.RTlog10 * pH

        # add the potential related to the ionic strength
        dG_prime -= 2.91482 * (z ** 2 - nH) * sqrt_I / (1 + 1.6 * sqrt_I)

        # add the potential related to the Mg ions
        if nMg > 0:
            dG_prime += nMg * \
                  (constants.RTlog10 * pMg - constants.MG_FORMATION_ENERGY)

        logging.debug('nH = %d, z = %d, dG0 = %.1f --> dG0\' = %.1f' %
                      (nH, z, dG, dG_prime))
        return dG_prime

    def __unicode__(self):
        return self.kegg_id

    def __str__(self):
        return "nH = %d, nMg = %d, z = %d, dG0_f = %.2f, phase = %s" %\
            (self.number_of_hydrogens, self.number_of_mgs, self.net_charge,
             self.formation_energy, self.phase)


class SpeciesGroup(models.Model):
    """
        A set of different species (AKA pseudoisomers/protonation states) that
        a compound can have. There is a separation between SpeciesGroup and
        Compound in order to associate a compound with more than one source
        of data, e.g. Alberty vs. Component Contribution
    """
    # The ID of the compound in KEGG.
    kegg_id = models.CharField(max_length=10)

    # The species in this group.
    species = models.ManyToManyField(Specie)

    # The priority of this group.
    priority = models.IntegerField()

    # The source of these values.
    formation_energy_source = models.ForeignKey(ValueSource,
                                                on_delete=models.CASCADE)

    def __init__(self, *args, **kwargs):
        super(SpeciesGroup, self).__init__(*args, **kwargs)
        self._all_species = None

    def __str__(self):
        s = "KEGG ID = %s\npriority = %d\nformation_energy_source = %s\n" % \
            (self.kegg_id, self.priority, self.formation_energy_source)
        return s + '\n'.join([str(sp) for sp in self.all_species])

    def GetSpecies(self):
        """Gets the list of Species, potentially caching."""
        if self._all_species is None:
            self._all_species = self.species.all()
        return self._all_species
    all_species = property(GetSpecies)

    def GetPhaseSpecies(self, phase=constants.DEFAULT_PHASE):
        """
            Get a list of all species corresponding to a certain phase
        """
        return [s for s in self.GetSpecies()
                if s.get_phase_display() == phase]

    def GetSpeciesWithoutMg(self):
        """Gets the list of species without Mg bound."""
        for s in self.all_species:
            if s.number_of_mgs == 0:
                yield s
    all_species_no_mg = property(GetSpeciesWithoutMg)

    def StashTransformedSpeciesEnergies(self, aq_params):
        """
            Stash the transformed species formation energy in each one.
        """
        for species in self.all_species:
            species.transformed_energy = species.Transform(aq_params)

    def DeltaG0Prime(self, aq_params,
                     phase=constants.DEFAULT_PHASE):
        """
            Get a deltaG estimate for this group of species.

            Args:
                pH: the pH to estimate at.
                pMg: the pMg to estimate at.
                ionic_strength: the ionic strength to estimate at.

            Returns:
                The estimated delta G in the given conditions or None.
        """
        species = self.GetPhaseSpecies(phase)

        if not species:
            logging.warning('No data for this compound (%s) in this '
                            'phase (%s) and priority (%d)' %
                            (self.kegg_id, phase, self.priority))
            return None

        logging.debug('Calculating dG0\' for %s' % (self.kegg_id))
        if phase == constants.DEFAULT_PHASE:
            # Compute per-species transforms, scaled down by R*T.
            scaled_transforms = [(-s.Transform(aq_params) / constants.RT)
                                 for s in species]

            # Numerical issues: taking a sum of exp(v) for |v| quite large.
            # Use the fact that we take a log later to offset all values by a
            # constant (the minimum value).
            total = scaled_transforms[0]
            for i in range(1, len(scaled_transforms)):
                total = numpy.logaddexp(total, scaled_transforms[i])
            dg0_prime = -constants.RT * total
        else:
            if len(species) > 1:
                logging.error('only aqueous phase can have multiple species')
                raise Http404
            else:
                # compounds which are not in solution do not need to be
                # transformed since they have only one specie and there is no
                # equilibrium between pseudoisomers in that phase
                dg0_prime = species[0].Transform(aq_params)

        logging.debug('total: dG0\' = %.1f' % dg0_prime)
        return dg0_prime

    def GetSourceReferenceLink(self):
        source_name = self.formation_energy_source
        try:
            source = ValueSource.objects.get(name=source_name)
            url = source.url
        except ValueSource.DoesNotExist:
            url = "/data_refs"
        url = url or "/data_refs"
        return url

    source_reference = property(GetSourceReferenceLink)


class Compound(models.Model):
    """
        A single compound.
    """
    # The ID of the compound in KEGG.
    kegg_id = models.CharField(max_length=10, null=True)

    # InChI representation of the compound.
    inchi = models.CharField(max_length=2048, null=True)

    # The main name of the compound.
    name = models.ForeignKey(CommonName, on_delete=models.CASCADE,
                             related_name='primary_name_of')

    # A list of common names of the compound, used for searching.
    common_names = models.ManyToManyField(CommonName)

    # If present, this name should always be used.
    preferred_name = models.CharField(max_length=250, null=True)

    # A remark about this compound.
    note = models.TextField(null=True)

    # A link for detailed remarks about this compound.
    details_link = models.URLField(null=True)

    # The chemical formula.
    formula = models.CharField(max_length=500, null=True)

    # The molecular mass.
    mass = models.FloatField(null=True)  # In Daltons.

    # The number of electrons.
    num_electrons = models.IntegerField(null=True)

    # The various estimates of Delta G for this compound.
    species_groups = models.ManyToManyField(SpeciesGroup)

    # The index of this compound in the Component Contribution
    # training data. This index is needed for calcuating the uncertainty.
    index = models.IntegerField(null=True)

    # Group vector encoded as a sparse vector in dictionary format.
    # Needed for calculating uncertainty.
    group_vector = models.TextField(null=True)

    # An explanation for when no DeltaG0 estimate is available.
    no_dg_explanation = models.CharField(max_length=2048,
                                         blank=True,
                                         null=True)

    # Thumbnail image for the structure of the compound
    thumbnail = models.TextField(blank=True)

    # A single static parser
    FORMULA_PARSER = formula_parser.FormulaParser()

    def __init__(self, *args, **kwargs):
        super(Compound, self).__init__(*args, **kwargs)
        self._all_species_groups = None
        self._species_group_to_use = None
        self._priority = None

    def GetSpeciesGroupPriorities(self):
        """
            Get the priorities of species groups available.
        """
        return [sg.priority for sg in self.all_species_groups]

    def SetSpeciesGroupPriority(self, priority):
        """
            Set the priority of the species group to use.
        """
        for sg in self.species_groups.all():
            if sg.priority == priority:
                self._species_group_to_use = sg
                logging.debug('Setting priority for %s to %d' %
                              (self.kegg_id, priority))
                break

    def SetHighestPriority(self):
        """
            Set the priority to the highest one.
        """
        ps = self.GetSpeciesGroupPriorities()
        if not ps:
            return

        self.SetSpeciesGroupPriority(min(ps))

    def HasData(self):
        """
            Has enough data to display.
        """
        return self.mass and self.formula

    def FirstName(self):
        """
            Return the 'first' name of this compound.

            If a 'preferred_name' is set, returns that. Otherwise, returns
            the first name in the list of common names. Presumes that the
            list of names is in some order.
        """
        if self.preferred_name:
            return self.preferred_name

        names = list(self.common_names.all())
        return names[0].name

    def NameSlug(self):
        """Return a name with no whitespace or dashes.

        Slug will also never begin with a number to avoid
        confusing the reaction parser.
        """
        slug = 'C_%s' % slugify(str(self.name))
        return slug.replace('-', '_')

    def DeltaG0Prime(self, aq_params,
                     phase=None):
        """
            Get a deltaG estimate for the given compound.

            Args:
                pH: the PH to estimate at.
                ionic_strength: the ionic strength to estimate at.
                temp: the temperature to estimate at.

            Returns:
                The estimated delta G in the given conditions or None.
        """
        phase = phase or self.GetDefaultPhaseName()
        sg = self._species_group
        if not sg:
            return None

        dg0_prime = sg.DeltaG0Prime(aq_params, phase)
        return dg0_prime

    def WriteStructureThumbnail(self):
        self.thumbnail = 'error'

        if self.inchi is None:
            return

        th_string = InChI2Thumbnail(str(self.inchi), output_format='png')
        if th_string is None:
            return

        self.thumbnail = base64.encodestring(th_string)

    def GetAtomBag(self):
        """
            Returns a dictionary of atoms and their counts for this compound.
        """
        if not self.formula:
            logging.debug('Formula is not defined for KEGG ID %s',
                          self.kegg_id)
            return None

        atom_bag = self.FORMULA_PARSER.GetAtomBag(self.formula)

        # Wildcards are not allowed.
        if 'R' in atom_bag:
            return None

        return atom_bag

    def GetNewPriorityLink(self, max_priority):
        return self.GetLink() + '&max_priority=%d' % max_priority

    def GetLink(self):
        """Returns the link to the stand-alone page for this compound."""
        if not self.kegg_id:
            return None
        return '/compound?compoundId=%s' % self.kegg_id

    def GetKeggLink(self):
        """Returns a link to the KEGG page for this compound."""
        if not self.kegg_id:
            return None

        return 'http://kegg.jp/dbget-bin/www_bget?cpd:%s' % self.kegg_id

    def GetSmallImageUrl(self):
        """Returns the URL of a small image of the compound structure."""
        if not self.kegg_id:
            return None

        return '/compound_image?compoundId=%s' % self.kegg_id

    def GetHtmlFormattedFormula(self):
        """Returns the chemical formula with HTML formatted subscripts."""
        if not self.formula:
            return None

        return re.sub(r'(\d+)', r'<sub>\1</sub>', self.formula)

    def GetSpeciesGroupToUse(self):
        """Gets the SpeciesGroup to use, potentially caching."""
        if not self._species_group_to_use:
            self.SetHighestPriority()
        return self._species_group_to_use

    def GetSpecies(self):
        """Gets the list of SpeciesFormationEnergies, potentially caching."""
        if self._species_group_to_use:
            return self._species_group_to_use.all_species

        # TODO(flamholz): Should we return something here?
        return None

    def SpeciesJson(self, species_group=None, aq_params=None,
                    phase=None):
        """Returns JSON for the species."""
        sg = species_group or self._species_group
        aq_params = aq_params or conditions.AqueousParams()
        phase = phase or self.GetDefaultPhaseName()
        dg0_prime = sg.DeltaG0Prime(aq_params, phase)
        logging.debug('dG0\' = %.1f ' % dg0_prime)
        d = {'source': str(sg.formation_energy_source),
             'dgzero_prime': {'value': round(dg0_prime, 1),
                              'pH': aq_params.pH,
                              'ionic_strength': aq_params.ionic_strength}}
        l = []
        for s in sg.all_species:
            l.append({'nh': int(s.number_of_hydrogens),
                      'charge': int(s.net_charge),
                      'nmg': int(s.number_of_mgs),
                      'dgzero': float(s.formation_energy)})
        d['species'] = l
        return d

    def AllSpeciesGroupsJson(self, aq_params, phase=constants.DEFAULT_PHASE):
        return [self.SpeciesJson(sg, aq_params, phase)
                for sg in self.all_species_groups]

    def ToJson(self, aq_params):
        d = {'name': str(self.name),
             'KEGG_ID': self.kegg_id,
             'InChI': self.inchi,
             'mass': self.mass,
             'formula': self.formula,
             'num_electrons': self.num_electrons,
             'thermodynamic_data': self.AllSpeciesGroupsJson(aq_params),
             'note': None}
        if self.note:
            d['note'] = self.note
        return d

    def ToCSVdG0Prime(self, priority=1, aq_params=None, phase=None):
        """
            returns a list of CSV rows with the following columns:
            kegg ID, name, dG0_prime, pH, ionic_strength, T, Note
        """
        phase = phase or self.GetDefaultPhaseName()
        aq_params = aq_params or conditions.AqueousParams()
        rows = []
        for sg in self.species_groups.filter(priority=priority):
            dG0_prime = sg.DeltaG0Prime(aq_params, phase=phase)
            if dG0_prime is not None:
                dG0_prime = dG0_prime.round(1)
            rows.append([self.kegg_id, self.name, dG0_prime,
                         aq_params.pH, aq_params.ionic_strength,
                         constants.DEFAULT_TEMP, None])
        return rows

    def ToCSVdG0(self, priority=1, phase=None):
        """
            returns a list of CSV rows with the following columns:
            kegg ID, name, dG0, nH, charge, nMg, Note
        """
        phase = phase or self.GetDefaultPhaseName()
        rows = []
        for sg in self.species_groups.filter(priority=priority):
            for s in sg.GetPhaseSpecies(phase):
                nH = int(s.number_of_hydrogens)
                charge = int(s.net_charge)
                nMg = int(s.number_of_mgs)
                dG0 = round(float(s.formation_energy), 1)
                rows.append([self.kegg_id, self.name, dG0, nH,
                             charge, nMg, None])
        return rows

    def GetSpeciesGroups(self):
        """Gets the list of SpeciesGroups."""
        if self._all_species_groups is None:
            self._all_species_groups = \
                self.species_groups.all().prefetch_related('species')

        return self._all_species_groups

    def HasSpeciesGroups(self):
        """Returns true if this compound has any species groups."""
        return self.species_groups.count() > 0

    def GetUniqueSpeciesGroups(self):
        """Iterator of unique species groups."""
        sources_set = set()
        for sg in self.all_species_groups:
            source = sg.formation_energy_source
            if source not in sources_set:
                sources_set.add(source)
                yield sg

    def GetPossiblePhaseNames(self):
        possible_phases = set()
        for sg in self.all_species_groups:
            for s in sg.all_species:
                possible_phases.add(s.phase)

        return sorted(possible_phases)

    def GetDefaultPhaseName(self):
        possible_phases = self.GetPossiblePhaseNames()
        if possible_phases == [] or constants.DEFAULT_PHASE in possible_phases:
            return constants.DEFAULT_PHASE
        else:
            return possible_phases[0]

    def GetDefaultPhase(self, conc=None):
        phase_name = self.GetDefaultPhaseName()
        return conditions._BaseConditions._GeneratePhase(phase_name, conc)

    def _GetDGSource(self):
        """Returns the source of the dG data."""
        if not self._species_group:
            return None

        return self._species_group.formation_energy_source

    def _GetSubstrateEnzymes(self):
        return self.substrate_for_enzymes.all()

    def _GetProductEnzymes(self):
        return self.product_of_enzymes.all()

    def _GetCofactorEnzymes(self):
        return self.cofactor_of_enzymes.all()

    def _GetSparseGroupVector(self):
        if self.group_vector is None:
            return None
        else:
            return json.loads(self.group_vector)

    _species_group = property(GetSpeciesGroupToUse)
    first_name = property(FirstName)
    name_slug = property(NameSlug)
    html_formula = property(GetHtmlFormattedFormula)
    link = property(GetLink)
    kegg_link = property(GetKeggLink)
    small_image_url = property(GetSmallImageUrl)
    all_common_names = property(lambda self: self.common_names.all())
    all_species = property(GetSpecies)
    all_species_groups = property(GetSpeciesGroups)
    has_species_groups = property(HasSpeciesGroups)
    unique_species_groups = property(GetUniqueSpeciesGroups)
    substrate_of = property(_GetSubstrateEnzymes)
    product_of = property(_GetProductEnzymes)
    cofactor_of = property(_GetCofactorEnzymes)
    dg0_prime = property(DeltaG0Prime)
    dg_source = property(_GetDGSource)
    sparse_gv = property(_GetSparseGroupVector)

    def StashTransformedSpeciesEnergies(self, aq_params):
        """Stash the transformed species formation energy in each one."""
        for sg in self.all_species_groups:
            sg.StashTransformedSpeciesEnergies(aq_params)

    def __unicode__(self):
        """Return a single string identifier of this Compound."""
        names = self.all_common_names
        if names:
            return str(names[0])
        return str(self.formula)

    @staticmethod
    def GetCompoundsByKeggId(kegg_ids):
        """Fetch compounds from a list of KEGG IDs.

        Args:
            kegg_ids: a list of KEGG IDs.

        Returns:
            A dictionary mapping KEGG ID to Compounds.
        """
        compounds = Compound.objects.select_related().filter(
            kegg_id__in=kegg_ids)
        return dict((c.kegg_id, c) for c in compounds if c is not None)


class Reactant(models.Model):
    """A compound and its coefficient."""
    # The compound.
    compound = models.ForeignKey(Compound, on_delete=models.CASCADE)

    # The coeff.
    coeff = models.FloatField(default=1.0)

    def __init__(self, *args, **kwargs):
        super(Reactant, self).__init__(*args, **kwargs)

    @staticmethod
    def GetOrCreate(kegg_id, coeff):
        """Attempts to fetch the CommonName object, or creates it if not present.

        Args:
            name: the name to use.

        Returns:
            A CommonName object.
        """
        try:
            r = Reactant.objects.get(compound__kegg_id=kegg_id,
                                     coeff=coeff)
            return r
        except Exception:
            c = Compound.objects.get(kegg_id=kegg_id)
            r = Reactant(compound=c, coeff=coeff)
            r.save()
            return r


class CompoundWithCoeff(models.Model):
    """A compound with a stoichiometric coefficient."""

    def __init__(self, coeff, compound, phase, name=None):
        """Construct a CompoundWithCoeff object.

        Args:
            coeff: the coefficient.
            compound: the Compound object.
            name: a string of the compound name.
            concentration: the concentration (molar).
        """
        self.compound = compound
        self.coeff = coeff
        self.phase = phase
        self._name = name

    def Clone(self):
        other = CompoundWithCoeff(self.coeff, self.compound,
                                  self.phase, self.name)
        return other

    def AbsoluteCoefficient(self):
        return numpy.abs(self.coeff)

    @staticmethod
    def FromReactant(reactant):
        return CompoundWithCoeff(reactant.coeff, reactant.compound,
                                 conditions.StandardAqueousPhase(),
                                 name=reactant.compound.FirstName())

    @staticmethod
    def FromId(coeff, kegg_id, phase=None, name=None):
        d = {'kegg_id': kegg_id, 'coeff': coeff, 'phase': phase}
        if name is not None:
            d['name'] = name
        return CompoundWithCoeff.FromDict(d)

    @staticmethod
    def FromDict(d):
        kegg_id = d['kegg_id']

        if 'compound' in d:
            compound = d['compound']
        else:
            try:
                compound = apps.get_model('gibbs.Compound').objects.get(kegg_id=kegg_id)
            except Exception:
                return None

        coeff = d.get('coeff', 1)
        name = d.get('name', compound.FirstName())
        phase_name = d.get('phase', None)
        if phase_name is None:
            phase_name = compound.GetDefaultPhaseName()
        conc = d.get('conc', None)
        phase = conditions._BaseConditions._GeneratePhase(phase_name, conc)

        logging.debug('Compound = %s, name = %s, coeff = %d, %s' %
                      (kegg_id, name, coeff, phase))
        return CompoundWithCoeff(coeff, compound, phase, name)

    def ResetConcentration(self):
        if not self.phase.IsConstant():
            self.phase.SetValue(None)

    def Minus(self):
        """Returns a new CompoundWithCoeff with coeff = -self.coeff."""
        return CompoundWithCoeff(-self.coeff, self.compound, self.phase,
                                 self.name)

    def __str__(self):
        name = self.name or str(self.compound)
        return '%d %s' % (self.coeff, name)

    def ToJson(self, include_species=True):
        d = {'coeff': self.coeff,
             'KEGG_ID': self.compound.kegg_id,
             'phase': str(self.phase),
             'name': str(self.compound.first_name),
             'source_used': None}

        if self.compound._species_group is not None:
            d['source_used'] = \
                str(self.compound._species_group.formation_energy_source)

        if include_species:
            d['species'] = self.compound.SpeciesJson()

        return d

    def _GetUrlParams(self):
        return ['reactantsId=%s' % self.kegg_id,
                'reactantsCoeff=%g' % self.coeff,
                'reactantsName=%s' % self.name,
                'reactantsPhase=%s' % self.phase.PhaseName(),
                'reactantsConcentration=%s' % self.phase.Value()]

    def GetCompoundList(self):
        sg = self.compound._species_group
        species = sg.GetPhaseSpecies(self.phase.PhaseName())
        logging.debug('KEGG ID = %s' % self.kegg_id)
        logging.debug('# species in phase %s = %d' %
                      (self.phase.PhaseName(), len(species)))

        l = []
        for s in species:
            l.append({'nh': int(s.number_of_hydrogens),
                      'charge': int(s.net_charge),
                      'nmg': int(s.number_of_mgs),
                      'dgzero': float(s.formation_energy)})
        logging.debug('compound_list = %s' % str(l))
        return l

    def GetName(self):
        """Gives a string name for this compound."""
        if self.compound.preferred_name:
            return self.compound.preferred_name
        if self._name:
            return self._name
        return str(self.compound.FirstName())

    def __eq__(self, other):
        """
            Check equality with another CompoundWithCoeff.

            Args:
                other: a second CompoundWithCoeff (or like object).
        """
        if self.coeff != other.coeff:
            return False

        if self.GetKeggID() != other.GetKeggID():
            return False

        if self.phase.Name() != other.phase.Name():
            return False

        return True

    def DeltaG0Prime(self, aq_params):
        dg0_prime = self.compound.DeltaG0Prime(aq_params,
                                               self.phase.PhaseName())
        if dg0_prime is None:
            return None
        return self.coeff * dg0_prime

    def GetPhaseName(self):
        return self.phase.PhaseName()

    def GetPossiblePhaseNames(self):
        return self.compound.GetPossiblePhaseNames()

    def HasMultiplePhases(self):
        return len(self.GetPossiblePhaseNames()) > 1

    def GetPhaseSubscript(self):
        # The compound C00288 represents the group of all carbonate species
        # and CO2 and is called CO2(total). It is technically aqueous, but
        # also contains a gas component, so the phase is not well defined
        # and hence the (total) subscript. In this case, we don't need to
        # add another subscript so we return an empty string
        if self.name == 'CO2(total)':
            return ''
        return self.phase.Subscript()

    def GetPhaseValueString(self):
        return '%.8g' % self.phase.HumanValueAndUnits()[0]

    def GetPhasePrefactor(self):
        return self.phase.HumanValueAndUnits()[1]

    def GetPhaseUnits(self):
        return self.phase.HumanValueAndUnits()[2]

    def GetPhaseIsConstant(self):
        return self.phase.IsConstant()

    def GetPhaseHumanString(self):
        value, prefactor, units = self.phase.HumanValueAndUnits_letters()
        return '%g %s' % (value, units)

    def GetKeggID(self):
        return self.compound.kegg_id

    name = property(GetName)
    abs_coeff = property(AbsoluteCoefficient)
    subscript = property(GetPhaseSubscript)
    phase_name = property(GetPhaseName)
    human_string = property(GetPhaseHumanString)
    phase_value_string = property(GetPhaseValueString)
    phase_prefactor = property(GetPhasePrefactor)
    phase_units = property(GetPhaseUnits)
    is_constant = property(GetPhaseIsConstant)
    possible_phases = property(GetPossiblePhaseNames)
    has_multiple_phases = property(HasMultiplePhases)
    kegg_id = property(GetKeggID)
