# -*- coding: utf-8 -*-
import logging
from gibbs import constants

class _BasePhase(object):
    def PhaseName(self):
        return NotImplementedError
    def Name(self):
        return NotImplementedError
    def Subscript(self):
        return NotImplementedError
    def Units(self):
        return NotImplementedError
    def Value(self):
        # the value in the given units (i.e. the number shown to the user)
        return 1
    def ValueString(self):
        return '%.2g' % self.Value()
    def IsConstant(self):
        return True
    def __str__(self):
        if self.Value() > 1e-2:
            return '%.2g %s' % (self.Value(), self.Units())
        
        if self.Value() > 1e-4:
            return '%.2g m%s' % (self.Value() * 1e3, self.Units())
        
        return '%2g μ%s' % (self.Value() * 1e6, self.Units())

class StandardAqueousPhase(_BasePhase):
    def PhaseName(self):
        return constants.AQUEOUS_PHASE_NAME
    def Name(self):
        return constants.STANDARD_AQUEOUS_PHASE_NAME
    def Subscript(self):
        return '(aq)'
    def Units(self):
        return 'M'

class StandardGasPhase(_BasePhase):
    def PhaseName(self):
        return constants.GAS_PHASE_NAME
    def Name(self):
        return constants.STANDARD_GAS_PHASE_NAME
    def Subscript(self):
        return '(g)'
    def Units(self):
        return 'bar'

class StandardLiquidPhase(_BasePhase):
    def PhaseName(self):
        return constants.LIQUID_PHASE_NAME
    def Name(self):
        return constants.STANDARD_LIQUID_PHASE_NAME
    def Subscript(self):
        return '(l)'
    def Units(self):
        return 'bar'

class StandardSolidPhase(_BasePhase):
    def PhaseName(self):
        return constants.SOLID_PHASE_NAME
    def Name(self):
        return constants.STANDARD_SOLID_PHASE_NAME
    def Subscript(self): 
        return '(s)'
    def Units(self):
        return 'bar'

class CustomAqueousPhase(StandardAqueousPhase):
    def __init__(self, concentration=1.0): # in units of M
        self._concentration = concentration
    def Name(self):
        return constants.CUSTOM_AQUEOUS_PHASE_NAME
    def IsConstant(self):
        return False
    def Value(self):
        return self._concentration

class CustomGasPhase(StandardGasPhase):
    def __init__(self, partial_pressure=1.0):
        self._partial_pressure = partial_pressure
    def Name(self):
        return constants.CUSTOM_GAS_PHASE_NAME
    def IsConstant(self):
        return False
    def Value(self):
        return self._partial_pressure
        
###############################################################################


class _BaseConditions(object):

    def __init__(self,
                 pH=constants.DEFAULT_PH,
                 pMg=constants.DEFAULT_PMG,
                 ionic_strength=constants.DEFAULT_IONIC_STRENGTH,
                 temperature=constants.DEFAULT_TEMP,
                 e_reduction_potential=constants.DEFAULT_ELECTRON_REDUCTION_POTENTIAL):
        self.pH = pH
        self.pMg = pMg
        self.ionic_strength = ionic_strength
        self.temperature = temperature
        self.e_reduction_potential = e_reduction_potential

    def __str__(self):
        raise NotImplementedError
    
    def _GetUrlParams(self):
        return ['conditions=%s' % self.__str__()]

    def GetTemplateDict(self):
        return {'ph': self.pH, 'pmg': self.pMg,
                'ionic_strength': self.ionic_strength,
                'temperature': self.temperature,
                'e_reduction_potential': self.e_reduction_potential,
                'conditions': self.__str__()}
    
    def GetPhase(self, kegg_id):
        raise NotImplementedError
        
class StandardConditions(_BaseConditions):
    
    def __str__(self):
        return constants.STANDARD_CONDITION_STRING
        
    def GetPhase(self, kegg_id):
        if kegg_id == 'C00001':
            return StandardLiquidPhase()
        else:
            return CustomAqueousPhase()

class MillimolarConditions(_BaseConditions):

    def __str__(self):
        return constants.MILLIMOLAR_CONDITION_STRING
    
    def GetPhase(self, kegg_id):
        if kegg_id == 'C00001':
            return StandardLiquidPhase()
        else:
            return CustomAqueousPhase(1e-3)
            
class CustomConditions(_BaseConditions):
    
    @staticmethod
    def _GetPhase(phase, value):
        if phase == constants.AQUEOUS_PHASE_NAME:
            return CustomAqueousPhase(value)
        if phase == constants.GAS_PHASE_NAME:
            return CustomGasPhase(value)
        if phase == constants.LIQUID_PHASE_NAME:
            return StandardLiquidPhase()
        if phase == constants.SOLID_PHASE_NAME:
            return StandardSolidPhase()    
        raise NotImplementedError
    
    def SetPhasesAndRatios(self, all_ids, all_phases, all_ratios):
        self._phases = {}
        for kegg_id, phase, ratio in zip(all_ids, all_phases, all_ratios):
            self._phases[kegg_id] = CustomConditions._GetPhase(phase, ratio)

    def __str__(self):
        return constants.CUSTOM_CONDITION_STRING
    
    def GetPhase(self, kegg_id):
        if kegg_id not in self._phases:
            logging.error('Condition requested for unknown id: %s', kegg_id)

        return self._phases[kegg_id]

    def _GetUrlParams(self):
        params = []
        params.append('conditions=%s' % self.__str__())
        for kegg_id, phase in self._phases.iteritems():
            params.append('reactantsPhase=%s' % phase.Name())
            params.append('reactantsConcentration=%s' % phase.Value())
        return params

###############################################################################

def GetConditions(name, all_ids=None, all_phases=None, all_ratios=None):
    
    if name == constants.STANDARD_CONDITION_STRING:
        return StandardConditions()

    if name == constants.MILLIMOLAR_CONDITION_STRING:
        return MillimolarConditions()
    
    if name == constants.CUSTOM_CONDITION_STRING:
        assert all_ids and all_phases and all_ratios
        cc = CustomConditions()        
        cc.SetPhasesAndRatios(all_ids, all_phases, all_ratios)
        return cc

    logging.error('unrecognized condition name: ' + name)
    return None