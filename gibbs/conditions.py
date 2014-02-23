import logging

class _BasePhase(object):
    def Name(self):
        return NotImplementedError
    def Subscript(self):
        return NotImplementedError
    def Units(self):
        return NotImplementedError
    def Value(self):
        return 1
    def ValueString(self):
        return '%.2g' % self.Value()
    def IsConstant(self):
        return True
    def GetRatio(self):
        return 1
    def __str__(self):
        return '%g %s' % (self.Value(), self.Units())

class StandardAqueousPhase(_BasePhase):
    def Name(self):
        return 'aqueous solution'
    def Subscript(self):
        return '(aq)'
    def Units(self):
        return 'M'

class StandardGasPhase(_BasePhase):
    def Name(self):
        return 'gas'
    def Subscript(self):
        return '(g)'
    def Units(self):
        return 'Pa'

class StandardLiquidPhase(_BasePhase):
    def Name(self):
        return 'liquid'
    def Subscript(self):
        return '(l)'
    def Units(self):
        return 'Pa'

class StandardSolidPhase(_BasePhase):
    def Name(self):
        return 'solid'
    def Subscript(self): 
        return '(s)'
    def Units(self):
        return 'Pa'

class MillimolarAqueousPhase(StandardAqueousPhase):
    def Value(self):
        return 1e-3
    def __str__(self):
        return '1 mM'

class CustomAqueousPhase(StandardAqueousPhase):
    def __init__(self, concentration=1.0): # in units of M
        self._concentration = concentration
    def IsConstant(self):
        return False
    def Value(self):
        return self._concentration

class CustomGasPhase(StandardGasPhase):
    def __init__(self, partial_pressure=1.0):
        self._partial_pressure = partial_pressure
    def IsConstant(self):
        return False
    def Value(self):
        return self._partial_pressure
        
###############################################################################

STANDARD_CONDITION_STRING = 'standard'
MILLIMOLAR_CONDITION_STRING = 'mM'
CUSTOM_CONDITION_STRING = 'custom'

class _BaseConditions(object):

    def __str__(self):
        raise NotImplementedError
    
    def _GetUrlParams(self):
        return ['conditions=%s' % self.__str__()]
    
    def GetPhase(self, kegg_id):
        raise NotImplementedError
        
class StandardConditions(_BaseConditions):
    
    def __str__(self):
        return STANDARD_CONDITION_STRING
        
    def GetPhase(self, kegg_id):
        if kegg_id == 'C00001':
            return StandardLiquidPhase()
        else:
            return CustomAqueousPhase()

class MillimolarConditions(_BaseConditions):

    def __str__(self):
        return MILLIMOLAR_CONDITION_STRING
    
    def GetPhase(self, kegg_id):
        if kegg_id == 'C00001':
            return StandardLiquidPhase()
        else:
            return MillimolarAqueousPhase()
            
class CustomConditions(_BaseConditions):
    
    def __str__(self):
        return CUSTOM_CONDITION_STRING

    def __init__(self, kegg_ids, concentrations):
        self._phases = {}
        for kegg_id, concentration in zip((kegg_ids, concentrations)):
            if kegg_id == 'C00001':
                self._phases[kegg_id] = StandardLiquidPhase()
            else:
                self._phases[kegg_id] = CustomAqueousPhase(concentration)
    
    def GetPhase(self, kegg_id):
        if kegg_id not in self._phases:
            logging.error('Condition requested for unknown id: %s',
                          kegg_id)

        return self._phases[kegg_id]

    def _GetUrlParams(self):
        params = []
        params.append('conditions=%s' % self.__str__())
        for kegg_id, phase in self._phases.iteritems():
            params.append('reactantId=%s' % kegg_id)
            params.append('phase=%s' % phase.Name())
            params.append('value=%s' % phase.Value())
        return params
        
def GetConditions(name, all_ids=None, all_concentrations=None):
    
    if name == STANDARD_CONDITION_STRING:
        return StandardConditions()

    if name == MILLIMOLAR_CONDITION_STRING:
        return MillimolarConditions()
    
    if name == CUSTOM_CONDITION_STRING:
        assert all_ids and all_concentrations
        return CustomConditions(all_ids, all_concentrations)

    logging.error('unrecognized condition name: ' + name)
    return None