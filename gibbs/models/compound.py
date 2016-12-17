#!/usr/bin/python

import hashlib
import json
import logging
import numpy
import re
import base64

from django.http import Http404
from django.db import models
from django.utils.text import slugify
from .. import constants
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
    formation_energy_source = models.ForeignKey(ValueSource, null=True)

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
    formation_energy_source = models.ForeignKey(ValueSource)

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

        def transform(x): x.Transform(aq_params)
        logging.debug('Calculating dG0\' for %s' % (self.kegg_id))
        if phase == constants.DEFAULT_PHASE:
            # Compute per-species transforms, scaled down by R*T.
            scaled_transforms = [(-transform(s) / constants.RT)
                                 for s in species]

            # Numerical issues: taking a sum of exp(v) for |v| quite large.
            # Use the fact that we take a log later to offset all values by a
            # constant (the minimum value).
            total = scaled_transforms[0]
            for i in xrange(1, len(scaled_transforms)):
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
                dg0_prime = transform(species[0])

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
    name = models.ForeignKey(CommonName, related_name='primary_name_of')

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
        """Return a name with no whitespace or dashes."""
        slug = slugify(unicode(self.name))
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
            return unicode(names[0])
        return unicode(self.formula)

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
    compound = models.ForeignKey(Compound)

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


class StoredReaction(models.Model):
    """A reaction stored in the database."""
    # The ID of this reaction in KEGG.
    kegg_id = models.CharField(max_length=10, null=True)

    # a JSON representation of reaction (i.e. the reactants and coefficients)
    reactants = models.TextField(null=True)

    # a hash string for fast lookup of enzyme names by reaction
    reaction_hash = models.CharField(max_length=128, db_index=True)

    # a cache for the ToString() function which takes too long due to DB comm
    reaction_string = models.TextField(null=True)

    # a cache for the Link() function which takes too long due to DB comm
    link = models.TextField(null=True)

    @staticmethod
    def FromJson(rd):
        return StoredReaction(kegg_id=rd['RID'],
                              reactants=json.dumps(rd['reaction']))

    def GetSparseRepresentation(self):
        sparse = {}
        for coeff, kegg_id in json.loads(self.reactants):
            sparse[kegg_id] = coeff
        if 'C00080' in sparse:  # ignore H+ in stored reactions
            sparse.pop('C00080')
        return sparse

    @staticmethod
    def _CompoundToString(kegg_id, coeff):
        try:
            compound = Compound.objects.get(kegg_id=kegg_id)
            name = compound.FirstName()
        except Exception as e:
            logging.warning('Cannot find the name for %s' % kegg_id)
            logging.warning(str(e))
            name = kegg_id

        if coeff == 1:
            return name
        else:
            return "%g %s" % (coeff, name)

    def ToString(self):
        """
            String representation.
        """
        # TODO: need to replace the KEGG IDs with the common names of the
        #       compounds
        left = []
        right = []
        for coeff, kegg_id in json.loads(self.reactants):
            if coeff < 0:
                left.append(StoredReaction._CompoundToString(kegg_id, -coeff))
            elif coeff > 0:
                right.append(StoredReaction._CompoundToString(kegg_id, coeff))
        return "%s = %s" % (' + '.join(left), ' + '.join(right))

    @staticmethod
    def HashableReactionString(sparse):
        """Return a hashable string for a biochemical reaction.

        The string fully identifies the biochemical reaction up to
        directionality. If it is equal to another reaction's string,
        then they have identical stoichiometry up to their directionality.

        Args:
            sparse: a dictionary whose keys are kegg_id and values are
                    stoichiometric coefficients
        """
        if len(sparse) == 0:
            return ''

        # sort according to KEGG ID and normalize the stoichiometric
        # coefficients such that the coeff of the reactant with the lowest
        # ID will be 1
        kegg_id_list = sorted(sparse.keys())
        if sparse[kegg_id_list[0]] == 0:
            raise Exception('One of the stoichiometric coefficients is 0')
        norm_factor = 1.0 / sparse[kegg_id_list[0]]
        return ' + '.join(['%g %s' % (norm_factor*sparse[kegg_id], kegg_id)
                           for kegg_id in kegg_id_list])

    @staticmethod
    def HashReaction(sparse):
        md5 = hashlib.md5()
        md5.update(StoredReaction.HashableReactionString(sparse))
        return md5.hexdigest()

    @staticmethod
    def GetAtpHydrolysisHash():
        atp_sparse = {'C00002': -1, 'C00001': -1, 'C00008': 1, 'C00009': 1}
        return StoredReaction.HashableReactionString(atp_sparse)

    @staticmethod
    def GetCO2HydrationHash():
        co2_sparse = {'C00011': -1, 'C00001': -1, 'C00288': 1}
        return StoredReaction.HashableReactionString(co2_sparse)

    def GetHashableReactionString(self):
        return StoredReaction.HashableReactionString(
            self.GetSparseRepresentation())

    def GetHash(self):
        return StoredReaction.HashReaction(self.GetSparseRepresentation())

    def GenerateHash(self):
        self.reaction_hash = self.GetHash()
        self.reaction_string = self.ToString()
        self.link = self.Link()

    def __str__(self):
        """String representation."""
        return self.ToString()

    def Link(self):
        """
            Returns a link to this reaction's page.
        """
        try:
            rxn = self.ToReaction()
            return rxn.GetHyperlink(self.ToString())
        except AttributeError:
            raise Exception('Cannot find one of the compounds in the database')

    def ToReaction(self, priority=1, aq_params=None):
        """
            returns a list of CSV rows with the following columns:
            kegg ID, dG0_prime, pH, ionic_strength, T, Note
        """
        from gibbs import reaction
        reactants = [reaction.CompoundWithCoeff.FromId(coeff, kegg_id)
                     for coeff, kegg_id in json.loads(self.reactants)]
        rxn = reaction.Reaction(reactants)
        return rxn


class Enzyme(models.Model):
    """A single enzyme."""
    # EC class enzyme.
    ec = models.CharField(max_length=10)

    # A list of common names of the enzyme, used for searching.
    common_names = models.ManyToManyField(CommonName)

    # List of reactions this enzyme catalyzes.
    reactions = models.ManyToManyField(StoredReaction)

    def HasData(self):
        """Checks if it has enough data to display."""
        return self.ec and self.reactions.all()

    def ToJson(self):
        """Returns as a JSON-friendly object."""
        return {'EC': self.ec,
                'name': str(self.first_name)}

    def Link(self):
        """Returns a link to this reactions page."""
        return '/enzyme?ec=%s' % self.ec

    def KeggLink(self):
        """Returns a link to the KEGG page for this enzyme."""
        return 'http://kegg.jp/dbget-bin/www_bget?ec:%s' % self.ec

    def BrendaLink(self):
        """Returns a link to the BRENDA page for this enzyme."""
        return 'http://www.brenda-enzymes.org/php/result_flat.php4?ecno=%s' \
            % self.ec

    def AllReactions(self):
        """Returns all the reactions."""
        return self.reactions.all()

    def FirstName(self):
        """The first name in the list of names."""
        return self.all_common_names[0]

    def __unicode__(self):
        """Return a single string identifier of this enzyme."""
        return unicode(self.FirstName())

    all_common_names = property(lambda self: self.common_names.all())
    all_cofactors = property(lambda self: self.cofactors.all())
    first_name = property(FirstName)
    all_reactions = property(AllReactions)
    kegg_link = property(KeggLink)
    brenda_link = property(BrendaLink)
    link = property(Link)
