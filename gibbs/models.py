#!/usr/bin/python

import hashlib
import logging
import numpy
import re
from scipy.misc import logsumexp

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
    """A common name of a compound."""
    name = models.CharField(max_length=500, db_index=True)
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
    """The source of a particular numeric value."""
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

        return dG
    
    def __unicode__(self):
        return self.kegg_id
    
    
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
            
    def GetSpecies(self):
        """Gets the list of Species, potentially caching."""
        if self._all_species is None:
            self._all_species = self.species.all()
        return self._all_species
    all_species = property(GetSpecies)
    
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

    def DeltaG(self, pH=constants.DEFAULT_PH,
               pMg=constants.DEFAULT_PMG,
               ionic_strength=constants.DEFAULT_IONIC_STRENGTH):
        """
            Get a deltaG estimate for this group of species.
            
            Args:
                pH: the pH to estimate at.
                pMg: the pMg to estimate at.
                ionic_strength: the ionic strength to estimate at.
            
            Returns:
                The estimated delta G in the given conditions or None.
        """
        if not self.all_species:
            # No data...
            return None
        
        # Compute per-species transforms, scaled down by R*T.
        transform = lambda x: x.Transform(pH=pH, pMg=pMg, ionic_strength=ionic_strength)
        scaled_transforms = [(-transform(s) / constants.RT)
                             for s in self.all_species]
        
        # Numerical issues: taking a sum of exp(v) for |v| quite large.
        # Use the fact that we take a log later to offset all values by a 
        # constant (the minimum value).
        return -constants.RT * logsumexp(scaled_transforms)
    

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
    
    # An explanation for when no DeltaG estimate is available.
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
        """Get the priorities of species groups available."""
        return [sg.priority for sg in self.all_species_groups]
    
    def SetSpeciesGroupPriority(self, priority):
        """Set the priority of the species group to use."""
        for sg in self.species_groups.all():
            if sg.priority == priority:
                self._species_group_to_use = sg
                break
    
    def SetHighestPriority(self):
        """Set the priority to the highest one."""
        ps = self.GetSpeciesGroupPriorities()
        if not ps:
            return
        
        self.SetSpeciesGroupPriority(min(ps))
    
    def HasData(self):
        """Has enough data to display."""
        return self.mass and self.formula
    
    def FirstName(self):
        """Return the 'first' name of this compound.
        
        If a 'preferred_name' is set, returns that. Otherwise, returns
        the first name in the list of common names. Presumes that the
        list of names is in some order.
        """
        if self.preferred_name:
            return self.preferred_name
        
        names = list(self.common_names.all())
        return names[0].name
    
    def DeltaG(self, pH=constants.DEFAULT_PH,
               pMg=constants.DEFAULT_PMG,
               ionic_strength=constants.DEFAULT_IONIC_STRENGTH):
        """Get a deltaG estimate for the given compound.
        
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
        
        return sg.DeltaG(pH, pMg, ionic_strength)

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
        """Returns a dictionary of atoms and their counts for this compound."""
        if not self.formula:
            logging.warning('Formula is not defined for KEGG ID %s', self.kegg_id)
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
        if self._species_group_to_use:
            return self._species_group_to_use
        
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
                    ionic_strength=constants.DEFAULT_IONIC_STRENGTH):
        """Returns JSON for the species."""
        sg = species_group or self._species_group
        l = []
        d = {'source': str(sg.formation_energy_source),
             'dgzero_tag': {'value': round(sg.DeltaG(pH, pMg, ionic_strength), 1),
                            'pH': pH,
                            'ionic_strength': ionic_strength}}
        for s in sg.all_species:
            l.append({'nh': int(s.number_of_hydrogens),
                      'charge': int(s.net_charge),
                      'nmg': int(s.number_of_mgs),
                      'dgzero': float(s.formation_energy)})
        d['species'] = l
        return d
    
    def AllSpeciesGroupsJson(self, pH=constants.DEFAULT_PH,
                             pMg=constants.DEFAULT_PMG,
                             ionic_strength=constants.DEFAULT_IONIC_STRENGTH):
        return [self.SpeciesJson(sg, pH=pH, pMg=pMg, ionic_strength=ionic_strength)
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
                      ionic_strength=constants.DEFAULT_IONIC_STRENGTH):
        """
            returns a list of CSV rows with the following columns:
            kegg ID, name, dG0_prime, pH, ionic_strength, T, Note
        """
        rows = []
        for sg in self.species_groups.filter(priority=priority):
            dG0_prime = round(sg.DeltaG(pH, pMg, ionic_strength), 1)
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
    dgf_zero_tag = property(DeltaG)
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
    
    # The list of substrates.
    substrates = models.ManyToManyField(Reactant, related_name='substrate_in')
    
    # The list of products.
    products = models.ManyToManyField(Reactant, related_name='product_in')

    # a hash string for fast lookup of enzyme names by reaction    
    hash = models.CharField(max_length=2048)
    
    @staticmethod
    def _SideString(side):
        """Returns a string representation for a single side of a reaction.
        
        Args:
            side: the list of CompoundWithCoeff objects representing the side.
        """
        l = []
        for r in side:
            if r.coeff == 1:
                l.append(r.compound.FirstName())
            else:
                l.append('%d %s' % (r.coeff,
                                    r.compound.FirstName()))
        return ' + '.join(l)

    def ReactionString(self):
        """Get the string representation of this reaction."""
        return '%s <=> %s' % (self._SideString(self.substrates.all()),
                              self._SideString(self.products.all()))
    
    def GetReactants(self):
        all_reactants = []
        for r in self.substrates.all():
            all_reactants.append((-r.coeff, r.compound.kegg_id, r.compound.FirstName()))
        for r in self.products.all():
            all_reactants.append((r.coeff, r.compound.kegg_id, r.compound.FirstName()))
        return all_reactants
    reactants = property(GetReactants)
    
    @staticmethod
    def HashableReactionString(substrates, products):
        """Return a hashable string for a biochemical reaction.
        
        The string fully identifies the biochemical reaction up to directionality.
        If it is equal to another reaction's string, then they have identical
        stoichiometry up to their directionality.
        
        Args:
            substrates: the substrates; a list of Reactants or like objects.
            products: the products; a list of Reactants or like objects.
        """
        sort_key = lambda r: r.compound.kegg_id
        make_str = lambda r: '%g %s' % (r.coeff, r.compound.kegg_id)
        is_not_hydrogen = lambda r: r.compound.kegg_id != 'C00080'
        
        substrates_strs = map(make_str,
                              sorted(filter(is_not_hydrogen, substrates),
                                     key=sort_key))
        rside_str = ' + '.join(substrates_strs)
        rside_hash = str(hash(rside_str))
        
        products_strs = map(make_str,
                            sorted(filter(is_not_hydrogen, products),
                                   key=sort_key))
        pside_str = ' + '.join(products_strs)
        pside_hash = str(hash(pside_str))
        
        sides = ['%s%s' % (rside_hash, rside_str),
                 '%s%s' % (pside_hash, pside_str)]
        sides.sort()
        return '%s <=> %s' % (sides[0], sides[1])
    
    def GetHashableReactionString(self):
        """Get a hashable string identifying this chemical reaction."""
        return self.HashableReactionString(self.substrates.all(),
                                           self.products.all())
    
    @staticmethod
    def HashReaction(substrates, products):
        md5 = hashlib.md5()
        md5.update(StoredReaction.HashableReactionString(substrates, products))
        return md5.hexdigest()
    
    def GetHash(self):
        """Returns a string hash of this reaction for easy identification."""
        return self.HashReaction(self.substrates.all(),
                                 self.products.all())
    
    def __hash__(self):
        """Makes stored reactions hashable."""
        return hash(self.GetHash())
    
    def __str__(self):
        """String representation."""
        return self.ReactionString()
    
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
        coeff_compound_pairs = []
        for reactant in self.substrates.all():
            coeff_compound_pairs.append((-reactant.coeff, reactant.compound))
        for reactant in self.products.all():
            coeff_compound_pairs.append((reactant.coeff, reactant.compound))
        
        dG0_prime = 0
        for coeff, compound in coeff_compound_pairs:
            dG0_f_prime = compound.DeltaG(pH, pMg, ionic_strength)
            if dG0_f_prime is not None:
                dG0_prime += coeff * dG0_f_prime
            else:
                return [(self.kegg_id, None, pH, ionic_strength,
                         constants.DEFAULT_TEMP,
                         "One of the compounds has no formation energy")]
        
        dG0_prime = round(dG0_prime, 1)
        return [(self.kegg_id, dG0_prime, pH, ionic_strength,
                 constants.DEFAULT_TEMP, None)]
        
    link = property(Link)
    reaction_string = property(ReactionString)

class ConservationLaw(models.Model):
    """ conservation laws which every reaction query must be checked against. """
    # The list of products.
    #reactants = models.ManyToManyField(Reactant, related_name='reactant_in')
    reactants = models.TextField(null=True)
    
    msg = models.TextField(null=True)
    
    def GetSparseRepresentation(self):
        sparse = {}
        for coeff, kegg_id in json.loads(self.reactants):
            sparse[kegg_id] = coeff
        return sparse
    
#    def __str__(self):
#        s_left = []
#        s_right = []
#        for reactant in self.reactants.all():
#            if abs(reactant.coeff) < 0.01:
#                continue
#            if reactant.coeff > 0:
#                if reactant.coeff == 1:
#                    s_right.append(reactant.compound.kegg_id)
#                else:
#                    s_right.append('%g %s' % (reactant.coeff, reactant.compound.kegg_id))
#            elif reactant.coeff < 0:
#                if reactant.coeff == -1:
#                    s_right.append(reactant.compound.kegg_id)
#                else:
#                    s_right.append('%g %s' % (-reactant.coeff, reactant.compound.kegg_id))
#        return ' + '.join(s_left) + ' => ' + ' + '.join(s_right)
    
class Enzyme(models.Model):
    """A single enzyme."""
    # EC class enzyme.
    ec = models.CharField(max_length=10)
    
    # A list of common names of the compound, used for searching.
    common_names = models.ManyToManyField(CommonName)
    
    # List of reactions this enzyme catalyzes.
    reactions = models.ManyToManyField(StoredReaction)
    
    # Compounds that this reaction
    substrates = models.ManyToManyField(Compound, related_name='substrate_for_enzymes')
    products = models.ManyToManyField(Compound, related_name='product_of_enzymes')
    cofactors = models.ManyToManyField(Compound, related_name='cofactor_of_enzymes')
    
    # TODO(flamholz): add more fields.
    
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
