from io import StringIO

from django.utils.text import slugify
from pathway.bounds import Bounds
from pathway.concs import ConcentrationConverter
from pathway.thermo_models import PathwayThermoModel

from pathway.max_min_driving_force import \
    PathwayParseError, InvalidReactionFormula, UnbalancedReaction, \
    ViolatesFirstLaw, ParsedPathway, ReactionMDFData, \
    CompoundMDFData, PathwayMDFData