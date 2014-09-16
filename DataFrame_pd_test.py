"""Unit tests for Dataframe_pd.py"""

import dataFrame_pd as df
import unittest
import numpy as np
import pandas as pd
import unitConversion as uc
import db_toolbox as SQL
import os

class Data:
    """This is the standard holder for all test data"""

    a = pd.Series([1.2,4.6,3.6,2.6], index = ['A', 'B', 'C', 'cheetah'])
    b = pd.Series([3.1,7.0,8.0,9.2], index = ['A', 'B', 'C', 'cheetah'])
    c = pd.Series([1.1,7.3,2.5,1.1], index = ['A', 'B', 'C', 'cheetah'])
    df1 = pd.DataFrame([a,b,c])


class CorrectInitialization(unittest.TestCase):
    """TESTS:
    +	Successfully create accurate dataframe from a dict of numpy arrays
    +   Successfully add units to that dict
    +   Raise error on non-existent column from unit dict
    +   Raise error if units are not given as strings
    """


    #This is a prototype right now
    def testCorrectInitialization(self):
        """The data in the data frame should match the data in the array_dict"""
        
        A = np.array([1.2,3.1,1.1])
        B = np.array([4.6,7.0,7.3])
        C = np.array([3.6,8.0,2.5])
        cheetah = np.array([2.6,9.2,1.1])

        data = df.Dataframe({'A':A,'B':B,'C':C,'cheetah':cheetah})
        
        #need a way to compare one dataframe to another
        self.assertTrue(np.all((data==Data.df1).values))

       

    def testUnitInitialization(self):
        """A unit dictionary must be appropriately applied to the dataframe"""
        units = {'A':'m/s', 'B':'kg/s', 'C':'Pa', 'cheetah':'s'}
        data = df.Dataframe([Data.a, Data.b, Data.c], units_dict = units) 

        self.assertEqual(units, data.units)

    def testNonexistentColumn(self):
        """Passing a column name in the units dict that isn't there should raise an error"""
        units = {'A':'m/s', 'B':'kg/s', 'C':'Pa', 'mouse':'s'}
        self.assertRaises(df.NoColumnError, df.Dataframe, [Data.a, Data.b, Data.c], units)

    def testBadUnitType(self):
        """Passing a non-string unit should raise an error"""
        units = {'A':'m/s', 'B':'kg/s', 'C':'Pa', 'cheetah':1.0}
        self.assertRaises(df.BadUnitError, df.Dataframe, [Data.a, Data.b, Data.c], units)
        
class ValueUnitTests(unittest.TestCase):
    """TESTS:
    +	Raise error if given column not in Dataframe
    +   Return a given value and unit together correctly
    +   Convert to a desired unit
    """
    def testNonexistentColumn(self):
        """Passing a column name that isn't there should raise an error"""
        units = {'A':'m/s', 'B':'kg/s', 'C':'Pa', 'cheetah':'s'}
        data = df.Dataframe([Data.a, Data.b, Data.c], units_dict = units)
        self.assertRaises(df.NoColumnError, data.val_units, 'mouse')

    def testCorrectOutput(self):
        """The dataframe should return a (value,unit) tuple correctly"""
        units = {'A':'m/s', 'B':'kg/s', 'C':'Pa', 'cheetah':'s'}
        data = df.Dataframe([Data.a, Data.b, Data.c], units_dict = units)
        vu = data.val_units('C')
        C = np.array([3.6,8.0,2.5])
        self.assertTrue((vu[0]==C).all())
        self.assertEqual(vu[1], 'Pa')

    def testUnitConversion(self):
        """The dataframe should return a (value,unit) tuple in the specified units"""
        units = {'A':'m/s', 'B':'kg/s', 'C':'Pa', 'cheetah':'s'}
        data = df.Dataframe([Data.a, Data.b, Data.c], units_dict = units)
        vu = data.val_units('A', 'ft/hr')
        A = np.array([1.2,3.1,1.1])
        conv = uc.UnitConverter()
        Ac = conv.convert_units(A, 'm/s', 'ft/hr')
        self.assertTrue((vu[0]==Ac).all())
        self.assertEqual(vu[1], 'ft/hr')

class LoadSQLTests(unittest.TestCase):
    """TESTS:
    + 	Raise error on bad interface
    +   Raise error on non-existent table
    +   Pass through query error
    +   Correctly load dataframe from database
    """

    def testBadInterface(self):
        """The dataframe should raise an error when a bad interface is passed"""
        data = df.Dataframe()
        self.assertRaises(df.dfSQLError, data.SQL_load_data, "dog")

    def testSQLError(self):
        """The dataframe should raise an error when the SQL library sends one back"""
        data = df.Dataframe()
        #error - no database selected
        interface = SQL.db_interface(host = "192.168.10.20", user = "chris", passwd = "cmp87ud01")
        interface.connect()
        self.assertRaises(df.dfSQLError, data.SQL_load_data, interface)

        q = SQL.use_Query("gas_unit_test")
        interface.query(q)

        #error - bad table
        self.assertRaises(df.dfSQLError, data.SQL_load_data, interface, "moat")        

        #error - bad condition
        self.assertRaises(df.dfSQLError, data.SQL_load_data, interface, "dataframe_pd_test_table", ["denver>2.0"])
        
    def testCorrectSQLLoad(self):
        """The dataframe should correctly be loaded from an SQL table"""
        data = df.Dataframe()
        interface = SQL.db_interface(host = "192.168.10.20", user = "chris", passwd = "cmp87ud01")
        interface.connect()
        q = SQL.use_Query("gas_unit_test")
        interface.query(q)

        data.SQL_load_data(interface, table = "dataframe_pd_test_table")
        self.assertTrue(np.all((data==Data.df1).values))

class UploadSQLTests(unittest.TestCase):
    """TESTS:
    +	Raise error on bad interface
    
    +   Sucessfuly load database w/new values
    +   Sucessully load database when dataframe contains nan
    +   Load updated data into SQL table successfully
    """
    def testBadInterface(self):
        """The dataframe should raise an error when a bad interface is passed"""
        
        data = df.Dataframe()
        self.assertRaises(df.dfSQLError, data.SQL_load_data, "dog")

    def testInsertLoad(self):
        """The dataframe should successfully add new records to an existing table"""
        # set up the database -- need to make a "safe" user for this, talk to Adrian about this
        os.system("mysql -h 192.168.10.20 -u chris -p < dataframe_pd_setup.sql")
        
        # go get the data from the database
        data = df.Dataframe()
        interface = SQL.db_interface(host = "192.168.10.20", user = "chris", passwd = "cmp87ud01")
        interface.connect()
        q = SQL.use_Query("gas_unit_test")
        interface.query(q)
        data.SQL_load_data(interface, table = "dataframe_pd_test_table")

        #upload the data to the other table
        data.SQL_upload_data(interface, table = "dataframe_upload_test_table")
        #download the data into a new dataframe
        data2 = df.Dataframe()
        data2.SQL_load_data(interface, table = "dataframe_upload_test_table")
         
        self.assertTrue(np.all((data==data2).values))

    def testInsertLoadNan(self):
        """The dataframe should successfully add new records to an existing table when NaN values are in the table"""
        # set up the database -- need to make a "safe" user for this, talk to Adrian about this
        os.system("mysql -h 192.168.10.20 -u chris -p < dataframe_pd_setup.sql")
        
        # go get the data from the database
        data = df.Dataframe()
        interface = SQL.db_interface(host = "192.168.10.20", user = "chris", passwd = "cmp87ud01")
        interface.connect()
        q = SQL.use_Query("gas_unit_test")
        interface.query(q)
        data.SQL_load_data(interface, table = "dataframe_pd_test_table")

        #replace a couple values with nan
        data.iloc[0][2] = np.nan
        data.iloc[1][1] = np.nan

        #upload the data to the other table
        data.SQL_upload_data(interface, table = "dataframe_upload_test_table")
        #download the data into a new dataframe
        data2 = df.Dataframe()
        data2.SQL_load_data(interface, table = "dataframe_upload_test_table")
         
        
        self.assertTrue(np.isnan(data2.iloc[0][2]))
        self.assertTrue(np.isnan(data2.iloc[1][1]))

    def testUpdateLoad(self):
        """The dataframe should successfully update existing records in an existing table"""
        # set up the database -- need to make a "safe" user for this, talk to Adrian about this
        os.system("mysql -h 192.168.10.20 -u chris -p < dataframe_pd_setup.sql")
        
        # go get the data from the database
        data = df.Dataframe()
        interface = SQL.db_interface(host = "192.168.10.20", user = "chris", passwd = "cmp87ud01")
        interface.connect()
        q = SQL.use_Query("gas_unit_test")
        interface.query(q)
        data.SQL_load_data(interface, table = "dataframe_pd_test_table")
        
        #create a new dataframe with one line
        
        data2 = df.Dataframe({'A':[1.2],'B':[1.0],'C':[2.0],'cheetah':[3.0]})
        
        data2.SQL_upload_data(interface, table = "dataframe_upload_test_table", index_col = 'A')
        
        #pull the data back into data
        data.SQL_load_data(interface, table = "dataframe_upload_test_table")

        #make sure the first row is equal
        self.assertEqual(1.2, data.iloc[0][0])
        self.assertEqual(1.0, data.iloc[0][1])
        self.assertEqual(2.0, data.iloc[0][2])
        self.assertEqual(3.0, data.iloc[0][3])

class CSVFileWriteTests(unittest.TestCase):
    """TESTS:
    +   Fail on bad filename
    0   Write a fully dataframe correctly
    0   Write a partial dataframe correctly
    +   Raise error when column not in dataframe
    0   Append to a CSV file correctly
    """
    def testBadFilename(self):
        """The dataframe should raise an error when a bad filename type is used"""
        A = np.array([1.2,3.1,1.1])
        B = np.array([4.6,7.0,7.3])
        C = np.array([3.6,8.0,2.5])
        cheetah = np.array([2.6,9.2,1.1])

        data = df.Dataframe({'A':A,'B':B,'C':C,'cheetah':cheetah})

        self.assertRaises(df.dfException, data.write_csv, 24.0)
        self.assertRaises(df.dfException, data.write_csv, None)

    def testCorrectWriteNewFull(self):
        """The dataframe should write to a new csv file correctly"""
        self.assertEqual(1,0)

    def testCorrectWritePartialFull(self):
        """The dataframe should write specified columsn to a new csv file correctly"""
        self.assertEqual(1,0)

    def testCorrectWriteAppend(self):
        """The dataframe should write to an existing csv file correctly"""
        self.assertEqual(1,0)

    def testNonexistentColumn(self):
        """Passing a column name that isn't there should raise an error"""
        A = np.array([1.2,3.1,1.1])
        B = np.array([4.6,7.0,7.3])
        C = np.array([3.6,8.0,2.5])
        cheetah = np.array([2.6,9.2,1.1])

        data = df.Dataframe({'A':A,'B':B,'C':C,'cheetah':cheetah})
        self.assertRaises(df.dfException, data.write_csv, "test.txt", "new", 'ninja')


class GlossaryReplaceTests(unittest.TestCase):
    """TESTS:
    +   Raise error on non-dict glossary
    +   Raise error on glossary column not in dataframe
    0   Raise error when value in glossary corresponds to an existing tag in the dataframe
    +   Correctly replace the names of the columns in the dataframe
    """

    def testBadGlossaryType(self):
        """Raise error on bad glossary type"""
        A = np.array([1.2,3.1,1.1])
        B = np.array([4.6,7.0,7.3])
        C = np.array([3.6,8.0,2.5])
        cheetah = np.array([2.6,9.2,1.1])

        data = df.Dataframe({'A':A,'B':B,'C':C,'cheetah':cheetah})
        self.assertRaises(df.BadGlossaryTypeError, data.glossary_replace, 2.0)
        self.assertRaises(df.BadGlossaryTypeError, data.glossary_replace, [2.3,4.2])
        self.assertRaises(df.BadGlossaryTypeError, data.glossary_replace, "fetid")

    def testNonexistentColumn(self):
        """Passing a column name that isn't there should raise an error"""
        A = np.array([1.2,3.1,1.1])
        B = np.array([4.6,7.0,7.3])
        C = np.array([3.6,8.0,2.5])
        cheetah = np.array([2.6,9.2,1.1])
        data = df.Dataframe({'A':A,'B':B,'C':C,'cheetah':cheetah})
        self.assertRaises(df.BadGlossaryTagError, data.glossary_replace, {'turtle':'BB'})

    def testDuplicateColumn(self):
        """Passing a column name as a replace value that already exists should raise an error"""
        A = np.array([1.2,3.1,1.1])
        B = np.array([4.6,7.0,7.3])
        C = np.array([3.6,8.0,2.5])
        cheetah = np.array([2.6,9.2,1.1])
        data = df.Dataframe({'A':A,'B':B,'C':C,'cheetah':cheetah})
        #self.assertRaises(df.BadGlossaryTagError, data.glossary_replace, {'cheetah':'A'})
        self.assertEqual(1,0)

    def testCorrectGlossaryReplace(self):
        """Correctly replace column names in the glossary"""
        A = np.array([1.2,3.1,1.1])
        B = np.array([4.6,7.0,7.3])
        C = np.array([3.6,8.0,2.5])
        cheetah = np.array([2.6,9.2,1.1])
        data = df.Dataframe({'A':A,'B':B,'C':C,'cheetah':cheetah})
        data.glossary_replace({'cheetah':'D'})
        data2 = df.Dataframe({'A':A, 'B':B,'C':C,'D':cheetah})
        self.assertTrue(np.all((data==data2).values))

class UnitSetTests(unittest.TestCase):
    """TESTS:
    +   Raise error on bad column name
    +   Raise error on bad unit name
    +   Correctly add units to the dataframe
    """
    def testNonexistentColumn(self):
        """Passing a column name that isn't there should raise an error"""
        A = np.array([1.2,3.1,1.1])
        B = np.array([4.6,7.0,7.3])
        C = np.array([3.6,8.0,2.5])
        cheetah = np.array([2.6,9.2,1.1])
        data = df.Dataframe({'A':A,'B':B,'C':C,'cheetah':cheetah})
        self.assertRaises(df.dfException, data.set_units,{'A':'m/s', 'frog':'m'})

    def testBadUnitType(self):
        """Passing a non-string unit should raise an error"""
        A = np.array([1.2,3.1,1.1])
        B = np.array([4.6,7.0,7.3])
        C = np.array([3.6,8.0,2.5])
        cheetah = np.array([2.6,9.2,1.1])
        data = df.Dataframe({'A':A,'B':B,'C':C,'cheetah':cheetah})
        self.assertRaises(df.dfException, data.set_units,{'A':'m/s', 'B':[2.3,5.3]})
	self.assertRaises(df.dfException, data.set_units,{'A':'m/s', 'B':1.2})

    def testCorrectUnitsAdd(self):
        """Unit dictionary should be added to the dataframe correctly"""
        A = np.array([1.2,3.1,1.1])
        B = np.array([4.6,7.0,7.3])
        C = np.array([3.6,8.0,2.5])
        cheetah = np.array([2.6,9.2,1.1])
        data = df.Dataframe({'A':A,'B':B,'C':C,'cheetah':cheetah})
        data.set_units({'A':'m/s', 'B':'m'})
        self.assertEqual(data.units, {'A':'m/s', 'B':'m'})


class FilterValsTests(unittest.TestCase):
    """TESTS:
    +   Raise error on bad types
    +   Raise error on bad actions
    +   Raise error on bad column
    +   Correctly filter high values
    +   Correctly filter low values
    """
    def testBadType(self):
        """Raise an error on bad unit types"""
        A = np.array([1.2,3.1,1.1])
        B = np.array([4.6,7.0,7.3])
        C = np.array([3.6,8.0,2.5])
        cheetah = np.array([2.6,9.2,1.1])
        data = df.Dataframe({'A':A,'B':B,'C':C,'cheetah':cheetah})
        self.assertRaises(df.BadArgumentError, data.filter_vals, 'cheetah', 'meep', 'high')
	self.assertRaises(df.BadArgumentError, data.filter_vals, 'cheetah', None, 'high')
	self.assertRaises(df.BadArgumentError, data.filter_vals, 'cheetah', [1.2, 3.4], 'high')

    def testBadAction(self):
        """Raise an error if the action is not in the action list"""
        A = np.array([1.2,3.1,1.1])
        B = np.array([4.6,7.0,7.3])
        C = np.array([3.6,8.0,2.5])
        cheetah = np.array([2.6,9.2,1.1])
        data = df.Dataframe({'A':A,'B':B,'C':C,'cheetah':cheetah})
        self.assertRaises(df.BadArgumentError, data.filter_vals, 'cheetah', 2.0, 'middle')

    def testNonexistentColumn(self):
        """Passing a column name that isn't there should raise an error"""
        A = np.array([1.2,3.1,1.1])
        B = np.array([4.6,7.0,7.3])
        C = np.array([3.6,8.0,2.5])
        cheetah = np.array([2.6,9.2,1.1])
        data = df.Dataframe({'A':A,'B':B,'C':C,'cheetah':cheetah})
        self.assertRaises(df.NoColumnError, data.filter_vals, 'frog', 2.0, 'high')

    def testCorrectFilterHigh(self):
        """Correctly filter values in high mode"""
        A = np.array([1.2,3.1,1.1])
        B = np.array([4.6,7.0,7.3])
        C = np.array([3.6,8.0,2.5])
        cheetah = np.array([2.6,9.2,1.1])
        data = df.Dataframe({'A':A,'B':B,'C':C,'cheetah':cheetah})
        filtered = data.filter_vals('cheetah', 2.0, 'high')
        
        self.assertTrue(np.isnan(filtered[0]))
        self.assertTrue(np.isnan(filtered[1]))
        self.assertEqual(filtered[2], 1.1)
        

    def testCorrectFilterLow(self):
        """Correctly filter values in low mode"""
        A = np.array([1.2,3.1,1.1])
        B = np.array([4.6,7.0,7.3])
        C = np.array([3.6,8.0,2.5])
        cheetah = np.array([2.6,9.2,1.1])
        data = df.Dataframe({'A':A,'B':B,'C':C,'cheetah':cheetah})
        filtered = data.filter_vals('A', 3.0, 'low')
        
        self.assertTrue(np.isnan(filtered[0]))
        self.assertTrue(np.isnan(filtered[2]))
        self.assertEqual(filtered[1], 3.1)

class ReplaceNans_w_NoneTests(unittest.TestCase):
    """TESTS:
    +   Raise error on bad column name
    +   Correctly replace nan values in the dataframe with None
    +   Do this for all columns, not just one
    """
    def testNonexistentColumn(self):
        """Passing a column name that isn't there should raise an error"""
        A = np.array([1.2,3.1,None])
        B = np.array([4.6,7.0,7.3])
        C = np.array([3.6,8.0,2.5])
        cheetah = np.array([2.6,9.2,1.1])
        data = df.Dataframe({'A':A,'B':B,'C':C,'cheetah':cheetah})
        self.assertRaises(df.NoColumnError, data.replace_None_w_nan, 'frogs')

    def testCorrectReplaceCol(self):
        """Correctly replace all NaN values in a dataframe column with Nones"""
        A = np.array([1.2,3.1,None])
        B = np.array([4.6,7.0,7.3])
        C = np.array([None,8.0,2.5])
        cheetah = np.array([2.6,9.2,1.1])
        data = df.Dataframe({'A':A,'B':B,'C':C,'cheetah':cheetah})
        data.replace_None_w_nan('A')
        self.assertTrue(np.isnan(data['A'][2]))
        self.assertEqual(None, data['C'][0])

    def testCorrectReplaceAll(self):
        """Correctly replace all NaN values in dataframe with Nones"""
        A = np.array([1.2,3.1,None])
        B = np.array([4.6,7.0,7.3])
        C = np.array([None,8.0,2.5])
        cheetah = np.array([2.6,9.2,1.1])
        data = df.Dataframe({'A':A,'B':B,'C':C,'cheetah':cheetah})
        data.replace_None_w_nan_all()
        self.assertTrue(np.isnan(data['A'][2]))
        self.assertTrue(np.isnan(data['C'][0]))
    

class FiniteValueTests(unittest.TestCase):
    """TESTS:
    +   Correctly return a column with only its finite values
    0   Correctly return a dataframe with only finite value columns
    +   Raise an error on a column not in the dataframe
    """

    def testCorrectFiniteValCol(self):
        """Correctly return a column of only finite values from a dataframe"""
        A = np.array([1.2,3.1,np.nan])
        B = np.array([4.6,7.0,7.3])
        C = np.array([np.nan,8.0,2.5])
        cheetah = np.array([2.6,9.2,1.1])
        data = df.Dataframe({'A':A,'B':B,'C':C,'cheetah':cheetah})
        f = data.finite_vals('A')
        
        self.assertEqual(len(f), 2)
        self.assertEqual(f[0], 1.2)
        self.assertEqual(f[1], 3.1)
        

    def testCorrectFiniteValDf(self):
        """Correctly return a dataframe of finite value columns for a given list of columns"""
        self.assertEqual(1,0)


    def testNonexistentColumn(self):
        """Passing a column name that isn't there should raise an error"""
        A = np.array([1.2,3.1,np.nan])
        B = np.array([4.6,7.0,7.3])
        C = np.array([np.nan,8.0,2.5])
        cheetah = np.array([2.6,9.2,1.1])
        data = df.Dataframe({'A':A,'B':B,'C':C,'cheetah':cheetah})
        self.assertRaises(df.NoColumnError, data.finite_vals, 'foo')
        


if __name__ == "__main__":
    unittest.main()
    
