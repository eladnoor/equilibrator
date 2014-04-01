#!/usr/bin/python

import hashlib
import logging
import numpy
import re
from scipy.misc import logsumexp

from django.http import Http404
from django.db import models
from gibbs import constants
from gibbs import formula_parser
#from django.core.files.base import ContentFile
import base64
import json

try:
    import indigo
    import indigo_renderer
    import openbabel
except ImportError:
    indigo = None
    indigo_renderer = None
    openbabel = None

class CommonName(models.Model):
    """
        A common name of a compound.
    """
    name = models.CharField(max_length=500)
    enabled = models.BooleanField(True)
    
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
    
    
class ValueSource(models.Model):
    """
        The source of a particular numeric value.
    """
    # The name of the source.
    name = models.CharField(max_length=100)
    
    # A citation name for the source. May be null.
    citation = models.CharField(max_length=4096, null=True)
    
    # The year of publication
    year = models.IntegerField()
    
    # The pubmed ID of the source if available.
    pubmed_id = models.CharField(max_length=128, null=True)

    # The DOI of the source if available.
    doi = models.CharField(max_length=128, null=True)
    
    # A link explaining the source.
    link = models.URLField(null=True)
    
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
        
    def Transform(self,
                  pH=constants.DEFAULT_PH,
                  pMg=constants.DEFAULT_PMG,
                  ionic_strength=constants.DEFAULT_IONIC_STRENGTH):
        """Transform this individual estimate to difference conditions."""

        dG = self.formation_energy
        
        # add the potential related to the pH
        if self.number_of_hydrogens > 0:
            dG += self.number_of_hydrogens * constants.RTlog10 * pH
        
        # add the potential related to the ionic strength
        dG -= 2.91482 * (self.net_charge ** 2 - self.number_of_hydrogens) * \
               numpy.sqrt(ionic_strength) / (1 + 1.6 * numpy.sqrt(ionic_strength))
        
        # add the potential related to the magnesium ions
        if self.number_of_mgs > 0:
            dG += self.number_of_mgs * \
                   (constants.RTlog10 * pMg - constants.MG_FORMATION_ENERGY)

        logging.debug('nH = %d, z = %d, dG0 = %.1f --> dG0\' = %.1f' % 
                      (self.number_of_hydrogens, self.net_charge,
                       self.formation_energy, dG))
        return dG
    
    def __unicode__(self):
        return self.kegg_id
    
    def __str__(self):
        return "nH = %d, nMg = %d, z = %d, dG0_f = %.2f, phase = %s" %\
            (self.number_of_hydrogens, self.number_of_mgs, self.net_charge,
             self.formation_energy, self.phase)
    
class SpeciesGroup(models.Model):
    """
        A set of different species (AKA pseudoisomers/protonation states) that
        a compound can have. There is a separation between SpeciesGroup and Compound
        in order to associate a compound with more than one source of data,
        e.g. Alberty vs. Component Contribution
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
        return s + '\n'.join([str(s) for s in self.all_species])

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

    def StashTransformedSpeciesEnergies(self, ph, pmg, ionic_strength):
        """
            Stash the transformed species formation energy in each one.
        """
        for species in self.all_species:
            species.transformed_energy = species.Transform(
                pH=ph, pMg=pmg, ionic_strength=ionic_strength)

    def DeltaG0Prime(self,
                     pH=constants.DEFAULT_PH,
                     pMg=constants.DEFAULT_PMG,
                     ionic_strength=constants.DEFAULT_IONIC_STRENGTH,
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
        
        transform = lambda x: x.Transform(pH=pH, pMg=pMg, 
                                          ionic_strength=ionic_strength)
        logging.debug('Calculating dG0\' for %s' % (self.kegg_id))
        if phase == constants.DEFAULT_PHASE:
            # Compute per-species transforms, scaled down by R*T.
            scaled_transforms = [(-transform(s) / constants.RT)
                                 for s in species]
            
            # Numerical issues: taking a sum of exp(v) for |v| quite large.
            # Use the fact that we take a log later to offset all values by a 
            # constant (the minimum value).
            if len(scaled_transforms) > 1:
                dg0_prime = -constants.RT * logsumexp(scaled_transforms)
            else:
                dg0_prime = -constants.RT * scaled_transforms[0]
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

    # Group vector encoded as a sparse vector in dictionary format.
    group_vector = models.TextField(null=True)
    
    # Replace this compound with another one.
    replace_with = models.ForeignKey('self', null=True)
    
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
    
    def DeltaG0Prime(self,
                     pH=constants.DEFAULT_PH,
                     pMg=constants.DEFAULT_PMG,
                     ionic_strength=constants.DEFAULT_IONIC_STRENGTH,
                     phase=constants.DEFAULT_PHASE):
        """
            Get a deltaG estimate for the given compound.
        
            Args:
                pH: the PH to estimate at.
                ionic_strength: the ionic strength to estimate at.
                temp: the temperature to estimate at.
            
            Returns:
                The estimated delta G in the given conditions or None.
        """
        sg = self._species_group
        if not sg:
            return None
        
        dg0_prime = sg.DeltaG0Prime(pH=pH, pMg=pMg, 
                                    ionic_strength=ionic_strength, 
                                    phase=phase)
        return dg0_prime

    def WriteStructureThumbnail(self, output_format='png'):
        if not indigo or not indigo_renderer or not openbabel:
            # the web server is not supposed to be running this function
            logging.error('Indigo/Openbabel are not installed, cannot draw structures.')

        if self.inchi is None:
            return

        _obConversion = openbabel.OBConversion()
        _obConversion.SetInFormat("inchi")
        _obConversion.SetOutFormat("smiles")
        obmol = openbabel.OBMol()
        _obConversion.ReadString(obmol, str(self.inchi))
        smiles = _obConversion.WriteString(obmol)

        _indigo = indigo.Indigo()
        _renderer = indigo_renderer.IndigoRenderer(_indigo)
        _indigo.setOption('render-output-format', output_format)
        _indigo.setOption('render-image-size', 250, 200)
        _indigo.setOption('render-margins', 10, 10)
        _indigo.setOption('render-stereo-style', 'none')
        _indigo.setOption('render-implicit-hydrogens-visible', False)
        _indigo.setOption('render-coloring', True)
        _indigo.setOption('render-bond-length', 50.0)
        _indigo.setOption('render-label-mode', 'hetero')

        try:
            indigo_mol = _indigo.loadMolecule(smiles)
            indigo_mol.aromatize()
            indigo_mol.layout()
            
            data = _renderer.renderToBuffer(indigo_mol).tostring()
            self.thumbnail = base64.encodestring(data)
        except indigo.IndigoException as e:
            logging.warning("Cannot draw structure of %s: %s" % (self.kegg_id,
                                                                         str(e)))

    def GetAtomBag(self):
        """
            Returns a dictionary of atoms and their counts for this compound.
        """
        if not self.formula:
            logging.debug('Formula is not defined for KEGG ID %s', self.kegg_id)
            return None
        
        atom_bag = self.FORMULA_PARSER.GetAtomBag(self.formula)
        
        # Wildcards are not allowed.
        if 'R' in atom_bag:
            return None
        
        return atom_bag

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
    
    def SpeciesJson(self, species_group=None,
                    pH=constants.DEFAULT_PH,
                    pMg=constants.DEFAULT_PMG,
                    ionic_strength=constants.DEFAULT_IONIC_STRENGTH,
                    phase=constants.DEFAULT_PHASE):
        """Returns JSON for the species."""
        sg = species_group or self._species_group
        dg0_prime = sg.DeltaG0Prime(pH=pH, pMg=pMg,
                                    ionic_strength=ionic_strength,
                                    phase=phase)
        logging.info(dg0_prime)
        d = {'source': str(sg.formation_energy_source),
             'dgzero_prime': {'value': round(dg0_prime, 1),
                              'pH': pH,
                              'ionic_strength': ionic_strength}}
        l = []
        for s in sg.all_species:
            l.append({'nh': int(s.number_of_hydrogens),
                      'charge': int(s.net_charge),
                      'nmg': int(s.number_of_mgs),
                      'dgzero': float(s.formation_energy)})
        d['species'] = l
        return d
    
    def AllSpeciesGroupsJson(self, pH=constants.DEFAULT_PH,
                             pMg=constants.DEFAULT_PMG,
                             ionic_strength=constants.DEFAULT_IONIC_STRENGTH,
                             phase=constants.DEFAULT_PHASE):
        return [self.SpeciesJson(sg, pH=pH, pMg=pMg, ionic_strength=ionic_strength,
                                 phase=phase)
                for sg in self.all_species_groups]
    
    def ToJson(self, pH=constants.DEFAULT_PH,
               pMg=constants.DEFAULT_PMG,
               ionic_strength=constants.DEFAULT_IONIC_STRENGTH):
        d = {'name': str(self.name), 
             'KEGG_ID': self.kegg_id,
             'InChI': self.inchi,
             'mass': self.mass,
             'formula': self.formula,
             'num_electrons': self.num_electrons,
             'thermodynamic_data': self.AllSpeciesGroupsJson(
                                pH=pH, pMg=pMg, ionic_strength=ionic_strength),
             'note': None}
        if self.note:
            d['note'] = self.note
        return d
    
    def ToCSVdG0Prime(self, priority=1,
                      pH=constants.DEFAULT_PH,
                      pMg=constants.DEFAULT_PMG,
                      ionic_strength=constants.DEFAULT_IONIC_STRENGTH,
                      phase=constants.DEFAULT_PHASE):
        """
            returns a list of CSV rows with the following columns:
            kegg ID, name, dG0_prime, pH, ionic_strength, T, Note
        """
        rows = []
        for sg in self.species_groups.filter(priority=priority):
            dG0_prime = round(sg.DeltaG0Prime(pH=pH, pMg=pMg, ionic_strength=ionic_strength, phase=phase), 1)
            rows.append([self.kegg_id, self.name, dG0_prime, pH, ionic_strength,
                         constants.DEFAULT_TEMP, None])
        return rows
            
    def ToCSVdG0(self, priority=1):
        """
            returns a list of CSV rows with the following columns:
            kegg ID, name, dG0, nH, charge, nMg, Note
        """
        rows = []
        for sg in self.species_groups.filter(priority=priority):
            for s in sg.all_species:
                nH = int(s.number_of_hydrogens)
                charge = int(s.net_charge)
                nMg = int(s.number_of_mgs)
                dG0 = round(float(s.formation_energy), 1)
                rows.append([self.kegg_id, self.name, dG0, nH, charge, nMg, None])
        return rows

    def GetSpeciesGroups(self):
        """Gets the list of SpeciesGroups."""
        if self._all_species_groups is None:
            self._all_species_groups = self.species_groups.all()
            
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
        if constants.DEFAULT_PHASE in possible_phases:
            return constants.DEFAULT_PHASE
        else:
            return possible_phases[0]

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
    
    _species_group = property(GetSpeciesGroupToUse)
    first_name = property(FirstName)
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
    dgf_zero_prime = property(DeltaG0Prime)
    dg_source = property(_GetDGSource)
    
    def StashTransformedSpeciesEnergies(self, ph, pmg, ionic_strength):
        """Stash the transformed species formation energy in each one."""
        for sg in self.all_species_groups:
            sg.StashTransformedSpeciesEnergies(ph, pmg, ionic_strength)
    
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
        return dict((c.kegg_id, c) for c in compounds if c != None)
    

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

    # a JSON representation of the reaction (i.e. the reactants and coefficients)
    reactants = models.TextField(null=True)
    
    # a hash string for fast lookup of enzyme names by reaction    
    reaction_hash = models.CharField(max_length=128, db_index=True)

    @staticmethod
    def FromJson(rd):
        return StoredReaction(kegg_id=rd['RID'],
                              reactants=json.dumps(rd['reaction']))
    
    def GetSparseRepresentation(self):
        sparse = {}
        for coeff, kegg_id in json.loads(self.reactants):
            sparse[kegg_id] = coeff
        if 'C00080' in sparse: # ignore H+ in stored reactions
            sparse.pop('C00080')
        return sparse
    
    @staticmethod
    def _CompoundToString(kegg_id, coeff):
        if coeff == 1:
            return kegg_id
        else:
            return "%g %s" % (coeff, kegg_id)

    def ToString(self):
        """
            String representation.
        """
        #TODO: need to replace the KEGG IDs with the common names of the compounds
        left = []
        right = []
        for coeff, kegg_id in json.loads(self.reactants):
            if coeff < 0:
                left.append(StoredReaction._CompoundToString(kegg_id, -coeff))
            elif coeff > 0:
                right.append(StoredReaction._CompoundToString(kegg_id, coeff))
        return "%s <=> %s" % (' + '.join(left), ' + '.join(right))
        
    @staticmethod
    def HashableReactionString(sparse):
        """Return a hashable string for a biochemical reaction.
        
        The string fully identifies the biochemical reaction up to directionality.
        If it is equal to another reaction's string, then they have identical
        stoichiometry up to their directionality.
        
        Args:
            sparse: a dictionary whose keys are kegg_id and values are 
                    stoichiometric coefficients
        """
        if len(sparse) == 0:
            return ''
        
        # sort according to KEGG ID and normalize the stoichiometric coefficients
        # such that the coeff of the reactant with the lowest ID will be 1
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
    
    def GetHashableReactionString(self):
        return StoredReaction.HashableReactionString(self.GetSparseRepresentation())
    
    def GetHash(self):
        return StoredReaction.HashReaction(self.GetSparseRepresentation())

    def GenerateHash(self):
        self.reaction_hash = self.GetHash()
    
    def __str__(self):
        """String representation."""
        return self.ToString()
    
    def Link(self):
        """Returns a link to this reaction's page."""
        return '/reaction?reactionId=%s' % self.kegg_id

    def ToCSVdG0Prime(self, priority=1,
                      pH=constants.DEFAULT_PH,
                      pMg=constants.DEFAULT_PMG,
                      ionic_strength=constants.DEFAULT_IONIC_STRENGTH):
        """
            returns a list of CSV rows with the following columns:
            kegg ID, dG0_prime, pH, ionic_strength, T, Note
        """
        dG0_prime = 0
        for kegg_id, coeff in json.loads(self.reactants):
            try:
                compound = models.Compound.objects.get(kegg_id=kegg_id)
                dG0_f_prime = compound.DeltaG0Prime(pH, pMg, ionic_strength)
                if dG0_f_prime is not None:
                    dG0_prime += coeff * dG0_f_prime
            except Exception:
                return [(self.kegg_id, None, pH, ionic_strength,
                        constants.DEFAULT_TEMP,
                        "One of the compounds has no formation energy")]
        
        return [(self.kegg_id, round(dG0_prime, 1), pH, ionic_strength,
                 constants.DEFAULT_TEMP, None)]
        
    link = property(Link)
    reaction_string = property(ToString)

class ConservationLaw(models.Model):
    """ conservation laws which every reaction query must be checked against. """
    # a JSON representation of the reaction (i.e. the reactants and coefficients)
    reactants = models.TextField(null=True)
    
    msg = models.TextField(null=True)
    
    def GetSparseRepresentation(self):
        sparse = {}
        for coeff, kegg_id in json.loads(self.reactants):
            sparse[kegg_id] = coeff
        return sparse
    
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
        return 'http://www.brenda-enzymes.org/php/result_flat.php4?ecno=%s' % self.ec

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
