#!/usr/bin/python

"""File for concentration maths."""


class ConcentrationConverter(object):

    # Units we might parse concentrations in.
    # Notice that the values are powers of 10 for the unit.
    UNITS_M  = 0
    UNITS_MM = -3
    UNITS_uM = -6

    UNITS_BY_NAME = {'MOLAR': UNITS_M,
                     'M': UNITS_M,
                     'MILLIMOLAR': UNITS_MM,
                     'MM': UNITS_MM,
                     'MICROMOLAR': UNITS_uM,
                     'UM': UNITS_uM}

    @classmethod
    def get_units(cls, unit_string):
        return cls.UNITS_BY_NAME.get(unit_string.upper())

    @classmethod
    def to_molar_units(cls, conc, from_units):
        factor = 10**from_units
        return conc * factor

    @classmethod
    def to_molar_string(cls, conc, from_units_string):
        from_units = cls.get_units(from_units_string)
        return cls.to_molar_units(conc, from_units)
