"""
tablib.dictionary.sbtab_dict
~~~~~~~~~~~~~
A wrapper object for handling an SBtab with multiple tables using
a dictinoary.
Also, includes methods for I/O between SQLite and SBtab.
"""
# -*- coding: utf-8 -*-
import tablib
from . import SBtabTools 
from . import tablibIO
from .SBtab import SBtabTable, SBtabError

class SBtabDict(dict):
    
    def __init__(self, sbtab_list):
        """
            Arguments:
                sbtab_list - a list of SBtabTable objects
        """
        self.fpath = ''
        self.sbtab_list = sbtab_list
        for m in sbtab_list:
            self[m.table_name] = m

    @staticmethod
    def FromSBtabFile(f):
        spreadsheet_file = tablibIO.haveTSV(f.read(), '\t')
        m = SBtabTools.oneOrMany(spreadsheet_file)
        sbtab_list = [SBtabTable(dset, 'dummy.tsv') for dset in m]
        sbtab_dict = SBtabDict(sbtab_list)
        sbtab_dict.fpath = 'dummy.tsv'
        return sbtab_dict

    @staticmethod
    def FromSBtab(fpath):
        spreadsheet_file = tablibIO.loadTSV(fpath, False)
        m = SBtabTools.oneOrMany(spreadsheet_file)
        sbtab_list = [SBtabTable(dset, fpath) for dset in m]
        sbtab_dict = SBtabDict(sbtab_list)
        sbtab_dict.fpath = fpath
        return sbtab_dict

    def GetColumnFromTable(self, table_name, column_name):
        """
            Returns:
                a list of the values in the column called 'column_name'
                in the table 'table_name'
        """
        column_index = self[table_name].columns_dict['!' + column_name]
        rows = self[table_name].getRows()
        return [r[column_index] for r in rows]

    def GetColumnsFromTable(self, table_name, column_names):
        """
            Arguments:
                table_name   - the name of the table in the SBtab file (without '!!')
                column_names - a list of column names from which to get the data (without '!')
                
            Returns:
                a list of lists containing the values corresponding to the
                columns in 'column_names' in the table 'table_name'
        """
        try:
            idxs = [self[table_name].columns_dict['!' + c] for c in column_names]
        except KeyError as e:
            all_columns = ', '.join(self[table_name].columns_dict.keys())
            raise KeyError('Cannot find the column %s in table "%s" in file %s. '
                           'Columns are: %s'
                           % (e, table_name, self.fpath, all_columns))
        return [map(r.__getitem__, idxs) for r in self[table_name].getRows()]
        
    def GetDictFromTable(self, table_name, key_column_name, value_column_name,
                         value_mapping=None):
        column_names = [key_column_name, value_column_name]
        keys, vals = zip(*self.GetColumnsFromTable(table_name, column_names))
        if value_mapping:
            return dict(zip(keys, map(value_mapping, vals)))
        else:
            return dict(zip(keys, vals))
        
    def GetTableAttribute(self, table_name, attribute_name):
        """
            Arguments:
                table_name     - the name of the table in the SBtab file (without '!!')
                attribute_name - a string with the attribute name
                
            Returns:
                A string containing the value of the attribute in that table,
                or None if the attribute does not exist
        """
        try:
            return self[table_name].getCustomTableInformation(attribute_name)
        except SBtabError:
            return None

