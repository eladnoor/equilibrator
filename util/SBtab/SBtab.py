"""
SBtab
=====

SBtab is a uniforming table format and designed for the use in
Systems Biology. Furthermore, it is a useful format to import stored information
into Python objects to process and manipulate it.

How to load tablib object:
Use "SBtabTools.openSBtab(filepath)" to create SBtab Python object.
Attention: Only 'tsv', 'csv', 'tab' and 'xls' are supported.

See specification for further information.
"""
#!/usr/bin/env python
import re
import copy
import tablib
from util.SBtab import tablibIO
from util.SBtab import misc

tables_without_name = []


class SBtabError(Exception):
    '''
    Base class for errors in the SBtab class.
    '''
    def __init__(self, message):
        self.message = message

    def __str__(self):
        return self.message


class SBtabTable():
    '''
    SBtabTable (version 0.9.0 06/10/2010)
    '''
    def __init__(self, table, filename):
        '''
        Creates SBtab Python object from tablib object.<

        Parameters
        ----------
        table : tablib object
            Containing one SBtab.
        filename : str
            Filename with extension.
        '''
        # Needed to be able to adress it from outside of the class for writing and reading
        self.filename = filename
        self.table = table

        #check if ascii stuff is violated
        if filename is not None:
            try:
                self.table = self.checkAscii(self.table)
            except:
                raise SBtabError('This is not a valid SBtab file. '
                                 'Try to check your file with the SBtab '
                                 'validator or read the SBtab specification.')

        # Delete tablib header to avoid complications
        if self.table.headers:
            self.table.headers = None

        # Create all necessary variables
        self.initializeTable()

        self.sbtab_dataset = []

    def initializeTable(self):
        '''
        Loads table informations and class variables.
        '''
        # Read the header row from table
        self.header_row = self._getHeaderRow()

        # Read the table information from header row
        (self.table_type, self.table_name, self.table_document,
         self.table_version, self.unique_key) = self.getTableInformation()
        
        # Read the columns of the table
        (self.columns, self.columns_dict, inserted_column, self.delimiter) = self.getColumns()

        # Read data rows
        self.value_rows = self.getRows(self.table_type, inserted_column)

        # Update the list and tablib object
        self.update()

    def checkAscii(self,table):
        '''
        Checks for ASCII violations, so that the parser will not crash if these occur.

        Parameters
        ----------
        table : str
            SBtab table as string.
        '''
        new_table = []

        for row in table:
            new_row = []
            for i,entry in enumerate(row):
                try:
                    print(entry)
                    new_row.append(str(entry).strip())
                except:
                    new_row.append('Ascii violation error! Please check input file!')
            new_table.append('\t'.join(new_row))

        tablibtable = tablibIO.importSetNew('\n'.join(new_table),
                                            self.filename + '.csv')

        return tablibtable

    def _getHeaderRow(self):
        '''
        Extracts the declaration row from the SBtab file.
        '''
        header_row = None
        # Find header row
        for row in self.table:
            for entry in row:
                if str(entry).startswith('!!'):
                    header_row = row
                    break
                elif str(entry).startswith('"!!'):
                    #rm1 = row.replace('""','#')
                    rm2 = row.remove('"')
                    header_row = rm2.replace('#','"')
                    break
        
        # Save string or raise error
        if not header_row:
            raise SBtabError('This is not a valid SBtab table, please use validator to check format or have a look in the specification!')
        else:
            header_row = ' '.join(header_row)

        # Replace double quotes by single quotes
        stupid_quotes = ['"','\xe2\x80\x9d','\xe2\x80\x98','\xe2\x80\x99','\xe2\x80\x9b','\xe2\x80\x9c','\xe2\x80\x9f','\xe2\x80\xb2','\xe2\x80\xb3','\xe2\x80\xb4','\xe2\x80\xb5','\xe2\x80\xb6','\xe2\x80\xb7']

        for squote in stupid_quotes:
            try: header_row = header_row.replace(squote, "'")
            except: pass
     
        # Split header row
        header_row = header_row.split(' ')

        # Delete spaces in header row
        while '' in header_row:
            header_row.remove('')

        header = ""
        for x in header_row[:-1]:
            header += x + ' '
        header += header_row[-1]

        return header

    def getTableInformation(self):
        '''
        Reads declaration row and stores the SBtab table attributes.
        '''

        # Initialize variables for unnamed table handling
        global tables_without_name
        no_name_counter = 0

        #header_row = self.getHeaderRow()
        # Save table type, otherwise raise error
        try: table_type = self.getCustomTableInformation('TableType')
        except: raise SBtabError('The TableType of the SBtab is not defined!')

        # Save table name, otherwise give name with number of unnamed tables
        try: table_name = self.getCustomTableInformation('TableName')
        #tn = re.search("TableName='([^']*)'", self.header_row)
        except:
            tables_without_name.append(table_type)
            for table_no_name in tables_without_name:
                if table_type == table_no_name:
                    no_name_counter = no_name_counter + 1
            table_name = table_type.capitalize() + '_' + str(no_name_counter)
            self.header_row += " TableName='" + table_name + "'"

        # Save table document, otherwise return None
        try: table_document = self.getCustomTableInformation('Document')
        except: table_document = None

        # save table version, otherwise return None
        try: table_version = self.getCustomTableInformation('SBtabVersion')
        except: table_version = None

        try: unique_key = self.getCustomTableInformation('UniqueKey')
        except: unique_key = 'True'

        return table_type, table_name, table_document, table_version, unique_key

    def getCustomTableInformation(self, attribute_name):
        '''
        Retrieves the value of a table attribute in the declaration line

        Parameters
        ----------
        attribute_name : str
           Name of the table attribute.
        '''
        if re.search("%s='([^']*)'" % attribute_name, self.header_row) != None:
            return re.search("%s='([^']*)'" % attribute_name, self.header_row).group(1)
        else:
            raise SBtabError('The %s of the SBtab is not defined!' % attribute_name)

    def getColumns(self):
        '''
        Extracts the column headers, adds mandatory first column name if necessary.
        '''
        # Save list of main columns
        for row in self.table:
            for entry in row:
                if str(row[0]).startswith('!') and not str(row[0]).startswith('!!'):
                    delimiter    = misc.getDelimiter(row)
                    column_names = list(row)
                    break

        # Insert mandatory first column if not existent
        inserted_column = False
        if not column_names[0].title() == '!' + self.table_type.title():
            column_names.insert(0, '!' + self.table_type.title())
            inserted_column = True

        # Get column positions
        columns = {}
        for i, column in enumerate(column_names):
            columns[column] = i

        return column_names, columns, inserted_column, delimiter

    def getRows(self, table_type='table', inserted=False):
        '''
        Extract the rows of the SBtab, add first column if necessary.

        Parameters
        ----------
        table_type : str (default 'table')
            Attribute TableType of the SBtab table. Only necessary, if first column was set automatically.
        inserted : Boolean
            True, if mandatory first column was set automatically.
        '''
        # Add row to list value_rows if row doesn't contain entries starting with '!'
        value_rows = []

        # Add to comments, if row starts with '%'
        self.comments = []

        for row in self.table:
            if str(row[0]).startswith('!'):            
                continue
            for i, entry in enumerate(row):
                #if str(entry).startswith('!'):
                #    break
                if str(entry).startswith('%'):
                    self.comments.append(list(row))
                    break
                else:
                    if len(row) == i + 1:
                        value_rows.append(list(row))

        # Insert value column if mandatory column was added automatically
        if inserted:
            if table_type == 'table':
                for i, row in enumerate(value_rows):
                    row.insert(0, 'TableRow_' + str(i + 1))
            else:
                for i, row in enumerate(value_rows):
                    row.insert(0, table_type[0].upper() + table_type[- 1].lower() + str(i + 1))

        return value_rows

    def changeValue(self, row, column, new):
        '''
        Change single value in the SBtab table by position in the table.

        Parameters
        ----------
        row : int
            Number of rows in the table. First row is number 1.
        column : int
            Number of columns in the table. First column is number 1.
        new : str
            New entry.
        '''
        self.value_rows[row - 1][column - 1] = new

        # Update object
        self.update()

    def changeValueByName(self, name, column_name, new):
        '''
        Change singe value in the SBtab by name of column and of the first row entry.

        Parameters
        ----------
        row : str
            Name of the entry in the first column.
        column_name : str
            Name of the column (without '!')
        new : str
            New entry.
        '''
        col = self.columns_dict['!' + column_name]
        for r in self.value_rows:
            if r[0] == name:
                r[col] = new

        # Update object
        self.update()

    def createList(self):
        '''
        Creates a list object of the SBtab Python object.
        '''
        # Create new list
        sbtab_list = []

        # Append the parts header row, main column row and value rows to the list
        sbtab_list.append(self.header_row)
        sbtab_list.append(self.columns)
        sbtab_list.append(self.value_rows)

        return sbtab_list

    def createDataset(self):
        '''
        Creates a tablib object of the SBtab Python object.
        '''
        # Initialize empty variables for conversion
        sbtab_temp = []
        self.sbtab_dataset = tablib.Dataset()

        # Create list of header
        header = [self.header_row]#.split(' ')

        # Delete spaces in header, main column and data rows
        header = [x.strip(' ') for x in header]
        self.columns = [x.strip(' ') for x in self.columns]
        for row in self.value_rows:
            try:
                for entry in row:
                    entry = entry.strip(' ')
            except:
                continue

        # Add header, main column and data rows to temporary list object
        sbtab_temp.append(header)
        sbtab_temp.append(self.columns)
        for row in self.value_rows:
            sbtab_temp.append(row)

        # Delete all empty entries at the end of the rows
        sb1 = []
        for row in sbtab_temp:
            if row[0] != '':
                sb1.append(row)
            #print '1. ',row
            #if len(row) > 1:
            #    while not row[-1]:
            #        del row[-1]

        # Make all rows the same length
        longest = max([len(x) for x in sb1])

        for row in sb1:
            if len(row) < longest:
                for i in range(longest - len(row)):
                    row.append('')
                self.sbtab_dataset.append(row)
            else:
                self.sbtab_dataset.append(row)

        return self.sbtab_dataset

    def addRow(self, row_list, position=None):
        '''
        Adds row to the table, if postion is None at the end of it.

        Parameters
        ----------
        row_list : list
            List of strings, containing the entries of the new row.
        position : int
            Position of new row in the table, 0 is on top.
        '''
        # Empty column to fill up sbtab_dataset with ''
        empty_list = []

        # Create temporary work copy
        sbtab_dataset = self.table

        # If new row is to small, add empty entries to new row
        if len(row_list) < len(sbtab_dataset.dict[0]):
            for i in range(len(sbtab_dataset.dict[0]) - len(row_list)):
                row_list.append('')

        # If new row is too long, add empty entries to sbtab_dataset
        elif len(row_list) > len(sbtab_dataset.dict[0]):
            for i in range(len(sbtab_dataset.dict[0])):
                empty_list.append('')
            for i in range(len(row_list) - len(sbtab_dataset.dict[0])):
                sbtab_dataset.rpush_col(empty_list)

        # If no position is set, add new row to the end
        if position is None:
            sbtab_dataset.rpush(row_list)
        else:
            sbtab_dataset.insert(position, row_list)

        # Update object
        self.table = sbtab_dataset
        self.initializeTable()

    def removeRow(self, position):
        '''
        Removes row from the table

        Parameters
        ----------
        position : int
            Position of row to be removed. Starting with 1.
        '''

        # Create temporary work copy
        sbtab_dataset = self.table

        del sbtab_dataset[position + 1]

        # Update object
        self.table = sbtab_dataset
        self.initializeTable()

    def addColumn(self, column_list, position=None):
        '''
        Adds a column to the table, if position is None at the end of it.

        Parameters
        ----------
        column_list : list
            List of strings, containing the entries of the new column.
        position : int
            Position of new column in the table, 0 is right.
        '''
        # Empty column to fill up sbtab_dataset with ''
        empty_list = []

        # If new column is to small, add empty entries to new column
        if len(column_list) < (len(self.sbtab_dataset.dict)-1):
            for i in range((len(self.sbtab_dataset.dict) - 1) - len(column_list)):
                column_list.append('')

        # If new column is to long, add empty entries to sbtab_dataset
        elif len(column_list) > (len(self.sbtab_dataset.dict) - 1):
            for i in range(len(self.sbtab_dataset.dict[0])):
                empty_list.append('')
            for i in range(len(column_list) - (len(self.sbtab_dataset.dict) - 1)):
                self.value_rows.append(empty_list)
                empty_list = copy.deepcopy(empty_list)

        # If no position is set, add new column to the end
        if not position:
            for i, row in enumerate(self.value_rows):
                row.append(column_list[i+1])
            self.columns_dict[column_list[0]] = len(self.columns)
            self.columns = self.columns_dict.keys()
        else:
            for i, row in enumerate(self.value_rows):
                row.insert(position - 1, column_list[i + 1])
            self.columns_dict[column_list[0]] = position - 1
            self.columns = self.columns_dict.keys()

        # Update object
        self.update()

    def removeColumn(self, position):
        '''
        Removes column from the table.

        Parameters
        ----------
        position : int
            Position of column to be removed. Sarting with 1.
        '''
        # Remove entries on position
        for row in self.value_rows:
            del row[position + 1]
        for column in self.columns_dict.keys():
            if self.columns_dict[column] == position - 1:
                del self.columns_dict[column]

        # Update object
        self.update()

    def writeSBtab(self, format_type, filename=None):
        '''
        Writes SBtab tablib object to file.

        Parameters
        ----------
        format_type : str
            File extension of the SBtab file. ('tsv', 'csv', 'tab', 'xls')
        filename : str
            Filename of the SBtab file without extension. Default is filename.
        '''
        if not filename:
            filename = self.filename[:-4]
        if format_type == 'tsv' or format_type == 'tab':
            tablibIO.writeTSV(self.sbtab_dataset, filename)
        elif format_type == 'csv':
            tablibIO.writeCSV(self.sbtab_dataset, filename)
        elif format_type == 'ods':
            tablibIO.writeODS(self.sbtab_dataset, filename)
        elif format_type == 'xls':
            tablibIO.writeXLS(self.sbtab_dataset, filename)
        else:
            raise SBtabError('The given file format is not supported: ' + format_type + '. Please use ".tsv", ".csv", ".tab" or ".xls" instead.')

    def duplicate(self):
        '''
        Creates a copy of the SBtab object.
        '''
        sbtab = copy.deepcopy(self)

        return sbtab

    def update(self):
        '''
        Updates the SBtab instance, list object, and tablib dataset.
        '''
        # Create tablib Dataset instance with new SBtab table
        self.table = self.createDataset()
        # Create list instance with new SBtab table
        self.sbtab_list = self.createList()

    def createSBtabDict(self):
        '''
        Creates a dict instance of the SBtab table.
        Keys are the column names, values are dicts. These contain the entries of the table.
        Keys are the entries in the first column, values are the current entries in the certain column.
        '''
        sbtab_dicts = {}
        for column_name in self.columns:
            sbtab_dicts[column_name] = {}
            for row in self.value_rows:
                sbtab_dicts[column_name][row[0]] = row[self.columns_dict[column_name]]

        return sbtab_dicts

    def transposeTable(self):
        '''
        Transposes SBtab table. Switches columns and rows.
        '''
        # Initialize new table data
        trans_columns = []
        trans_columns_dict = {}
        trans_value_rows = []

        # Save old table data
        columns = self.columns
        value_rows = self.value_rows

        # Append first entry to new column
        trans_columns.append(columns.pop(0))

        # Set new rows
        for column in columns:
                trans_value_rows.append([column])

        # Set new values in tables
        for row in value_rows:
            trans_columns.append(row.pop(0))
            for i, entry in enumerate(row):
                trans_value_rows[i].append(entry)

        # Write new columns dict
        for i, column in enumerate(trans_columns):
            trans_columns_dict[column] = i

        # Overwrite old table data
        self.columns = trans_columns
        self.columns_dict = trans_columns_dict
        self.value_rows = trans_value_rows

        self.update()

    def toDataFrame(self):
        import pandas as pd
        column_names = list(map(lambda s: s[1:], self.columns))
        df = pd.DataFrame(data=self.getRows(), columns=column_names)
        return df
    
    @staticmethod
    def fromDataFrame(df, document_name, table_type, table_name,
                      document, unit, sbtab_version='1.0'):
        table = tablib.Dataset()
        header = "!!SBtab DocumentName='%s' TableType='%s' TableName='%s' Document='%s' Unit='%s' SBtabVersion='%s'" % \
            (document_name, table_type, table_name, document, unit, sbtab_version)
        table.rpush([header] + [''] * (df.shape[1]-1))
        table.rpush(list(map(lambda s: '!' + s, df.columns)))
        for i, row in df.iterrows():
            table.append(row.tolist())
        return SBtabTable(table, None)
        