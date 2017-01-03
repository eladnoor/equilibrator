import logging
import re
import pyparsing

ELEMENTS = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na',
            'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti',
            'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge',
            'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo',
            'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te',
            'I', 'Xe', 'Cs', 'Ba', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os',
            'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn',
            'Fr', 'Ra', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds',
            'Rg', 'Uub', 'Uut', 'Uuq', 'Uup', 'Uuh', 'Uus', 'Uuo', 'La',
            'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho',
            'Er', 'Tm', 'Yb', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am',
            'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No',
            'R']  # R is not an alkane group, not an element but we allow it...


class _ElementAndCoeff(object):

    def __init__(self, el, coeff):
        self.el = el
        self.coeff = coeff

    def IsBlock(self):
        return False


class _Block(object):

    def __init__(self, elements):
        self.elements = elements

    def IsBlock(self):
        return True


class FormulaParser(object):

    def __init__(self):
        element = pyparsing.oneOf(ELEMENTS)
        coeff = pyparsing.Word(pyparsing.nums)
        optional_coeff = pyparsing.Optional(coeff)

        element_and_count = pyparsing.Forward()
        element_and_count << (element + optional_coeff)
        element_and_count.setParseAction(self.HandleElementAndCount)

        n_block = pyparsing.Forward()
        n_block << ('(' + pyparsing.OneOrMore(element_and_count) + ')n')
        n_block.setParseAction(self.HandleNBlock)

        element_or_block = pyparsing.Or([element_and_count, n_block])
        self.formula_parser = pyparsing.OneOrMore(element_or_block)

        self.formula_re = re.compile(r'\(?([A-Z][a-z]?)([0-9]*)(\)n)?')

    @staticmethod
    def HandleElementAndCount(block):
        if len(block) == 1:
            return _ElementAndCoeff(block[0], 1)

        element, count = block
        return _ElementAndCoeff(element, int(count))

    @staticmethod
    def HandleNBlock(block):
        return _Block(block[1:-1])

    def GetAtomBag(self, formula):
        """Returns a dictionary mapping atoms to counts for the formula."""
        if not formula:
            logging.error('Invalid formula %s', formula)
            return None

        atom_bag = {}
        results = self.formula_parser.parseString(formula)
        for item in results:
            if item.IsBlock():
                for el in item.elements:
                    coeff = 100*el.coeff
                    atom_bag[el.el] = atom_bag.setdefault(el.el, 0) + coeff
            else:
                atom_bag[item.el] = atom_bag.setdefault(item.el, 0) + \
                    item.coeff

        return atom_bag
