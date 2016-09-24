#!/usr/bin/python

"""File for concentration maths."""


class NoSuchUnits(Exception):
    pass


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
        """Performs the math of unit conversion.

        Args:
            conc: the concentration in "from_units" units.
            from_units: the source units, defined as the
                order-of-magnitude difference from molar
                units (see class definition above).
        
        Returns:
            The concentration in Molar.
        """
        factor = 10**from_units
        return conc * factor

    @classmethod
    def to_molar_string(cls, conc, from_units_string):
        """Performs the math of unit conversion.

        Fails with an unspecified error if from_unit_string
        is not specifying a recognized unit.
        
        Args:
            conc: the concentration in "from_units" units.
            from_units: the source units as a string.
        
        Returns:
            The concentration in Molar.
        """
        # TODO: handle case of unrecognized units.
        from_units = cls.get_units(from_units_string)
        if from_units is None:
            raise NoSuchUnits('"%s" could not be parsed' % from_units_string)
        return cls.to_molar_units(conc, from_units)
