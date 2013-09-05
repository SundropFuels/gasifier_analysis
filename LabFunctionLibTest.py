"""LabFunctionLibTest.py
   Unit tests for objects within the LabFunctionLib.py library
   Chris Perkins
   2012-01-24

   Version history:
   Version 0.1 - Tests for first porting over from R libraries

"""

import unittest
import LabFunctionLib as lfl
import numpy as np
import datetime
import db_toolbox as db
import unitConversion as uc
import dataFrame_v2 as df
import Cantera as ct

class LoadDataTests(unittest.TestCase):
    timestamp = np.ndarray(0, dtype = 'object')
    timestamp = np.append(timestamp, datetime.datetime(1981,7,6,13,13,12))
    timestamp = np.append(timestamp, datetime.datetime(1981,7,6,13,13,13))
    timestamp = np.append(timestamp, datetime.datetime(1981,7,6,13,13,14))
    timestamp = np.append(timestamp, datetime.datetime(1981,7,6,13,13,15))
    timestamp = np.append(timestamp, datetime.datetime(1981,7,6,13,13,16))

    ME_101 = np.ndarray(0, dtype = 'float64')
    ME_101 = np.append(ME_101, 3.5)
    ME_101 = np.append(ME_101, 3.4)
    ME_101 = np.append(ME_101, 3.6)
    ME_101 = np.append(ME_101, 3.5)
    ME_101 = np.append(ME_101, 3.5)

    TE_101 = np.ndarray(0, dtype = 'float64')
    TE_101 = np.append(TE_101, 25.0)
    TE_101 = np.append(TE_101, 26.0)
    TE_101 = np.append(TE_101, 25.0)
    TE_101 = np.append(TE_101, 27.0)
    TE_101 = np.append(TE_101, 25.0)

    TE_102 = np.ndarray(0, dtype = 'float64')
    TE_102 = np.append(TE_102, 120)
    TE_102 = np.append(TE_102, 121)
    TE_102 = np.append(TE_102, 122)
    TE_102 = np.append(TE_102, 121)
    TE_102 = np.append(TE_102, 122)

    PT_101 = np.ndarray(0, dtype = 'float64')
    PT_101 = np.append(PT_101, 2.51)
    PT_101 = np.append(PT_101, 2.50)
    PT_101 = np.append(PT_101, 2.55)
    PT_101 = np.append(PT_101, 2.46)
    PT_101 = np.append(PT_101, 2.49)

    PT_102 = np.ndarray(0, dtype = 'float64')
    PT_102 = np.append(PT_102, 1.01)
    PT_102 = np.append(PT_102, 1.11)
    PT_102 = np.append(PT_102, 1.21)
    PT_102 = np.append(PT_102, 1.00)
    PT_102 = np.append(PT_102, 1.09)

    PT_103 = np.ndarray(0, dtype = 'float64')
    PT_103 = np.append(PT_103, 0.33)
    PT_103 = np.append(PT_103, 0.34)
    PT_103 = np.append(PT_103, 0.32)
    PT_103 = np.append(PT_103, 0.33)
    PT_103 = np.append(PT_103, 0.34)

    MFC_101 = np.ndarray(0, dtype = 'float64')
    MFC_101 = np.append(MFC_101, 0.0)
    MFC_101 = np.append(MFC_101, 0.0)
    MFC_101 = np.append(MFC_101, 0.0)
    MFC_101 = np.append(MFC_101, 0.0)
    MFC_101 = np.append(MFC_101, 0.0)

    MFC_102 = np.ndarray(0, dtype = 'float64')
    MFC_102 = np.append(MFC_102, 7.1)
    MFC_102 = np.append(MFC_102, 7.09)
    MFC_102 = np.append(MFC_102, 7.11)
    MFC_102 = np.append(MFC_102, 7.10)
    MFC_102 = np.append(MFC_102, 7.10)

    MFC_103 = np.ndarray(0, dtype = 'float64')
    MFC_103 = np.append(MFC_103, 0.51)
    MFC_103 = np.append(MFC_103, 0.50)
    MFC_103 = np.append(MFC_103, 0.50)
    MFC_103 = np.append(MFC_103, 0.49)
    MFC_103 = np.append(MFC_103, 0.50)

    MFC_104 = np.ndarray(0, dtype = 'float64')
    MFC_104 = np.append(MFC_104, 1.01)
    MFC_104 = np.append(MFC_104, 1.00)
    MFC_104 = np.append(MFC_104, 0.97)
    MFC_104 = np.append(MFC_104, 0.99)
    MFC_104 = np.append(MFC_104, 1.05)

    FE_101 = np.ndarray(0, dtype = 'float64')
    FE_101 = np.append(FE_101, 51.0)
    FE_101 = np.append(FE_101, 51.3)
    FE_101 = np.append(FE_101, 50.9)
    FE_101 = np.append(FE_101, 50.8)
    FE_101 = np.append(FE_101, 51.2)
    
    CO_NDIR = np.ndarray(0, dtype = 'float64')
    CO_NDIR = np.append(CO_NDIR, 0.26)
    CO_NDIR = np.append(CO_NDIR, 0.27)
    CO_NDIR = np.append(CO_NDIR, 0.26)
    CO_NDIR = np.append(CO_NDIR, 0.27)
    CO_NDIR = np.append(CO_NDIR, 0.28)

    CO2_NDIR = np.ndarray(0, dtype = 'float64')
    CO2_NDIR = np.append(CO2_NDIR, 0.09)
    CO2_NDIR = np.append(CO2_NDIR, 0.089)
    CO2_NDIR = np.append(CO2_NDIR, 0.088)
    CO2_NDIR = np.append(CO2_NDIR, 0.089)
    CO2_NDIR = np.append(CO2_NDIR, 0.091)

    CH4_NDIR = np.ndarray(0, dtype = 'float64')
    CH4_NDIR = np.append(CH4_NDIR, 0.042)
    CH4_NDIR = np.append(CH4_NDIR, 0.041)
    CH4_NDIR = np.append(CH4_NDIR, 0.042)
    CH4_NDIR = np.append(CH4_NDIR, 0.041)
    CH4_NDIR = np.append(CH4_NDIR, 0.040)

    H2_GC = np.ndarray(5, dtype = 'float64')
    H2_GC[0] = 0.41
    H2_GC[1] = 0.40
    H2_GC[2] = 0.42
    H2_GC[3] = 0.42
    H2_GC[4] = 0.42

    CO_GC = np.ndarray(5, dtype = 'float64')
    CO_GC[0] = 0.26
    CO_GC[1] = 0.27
    CO_GC[2] = 0.26
    CO_GC[3] = 0.27
    CO_GC[4] = 0.28

    CO2_GC = np.ndarray(5, dtype = 'float64')
    CO2_GC[0] = 0.09
    CO2_GC[1] = 0.089
    CO2_GC[2] = 0.088
    CO2_GC[3] = 0.089
    CO2_GC[4] = 0.091

    CH4_GC = np.ndarray(5, dtype = 'float64')
    CH4_GC[0] = 0.042
    CH4_GC[1] = 0.041
    CH4_GC[2] = 0.042
    CH4_GC[3] = 0.041
    CH4_GC[4] = 0.040

    C2H6_GC = np.ndarray(5, dtype = 'float64')
    C2H6_GC[0] = 0.010
    C2H6_GC[1] = 0.010
    C2H6_GC[2] = 0.011
    C2H6_GC[3] = 0.010
    C2H6_GC[4] = 0.011

    N2_GC = np.ndarray(5, dtype = 'float64')
    N2_GC[0] = 0.186
    N2_GC[1] = 0.179
    N2_GC[2] = 0.168
    N2_GC[3] = 0.160
    N2_GC[4] = 0.1487

    Ar_GC = np.ndarray(5, dtype = 'float64')
    Ar_GC[0] = 0.0164
    Ar_GC[1] = 0.0112
    Ar_GC[2] = 0.0105
    Ar_GC[3] = 0.010
    Ar_GC[4] = 0.00929

    Counter = np.ndarray(5, dtype = 'int64')
    Counter[0] = 0
    Counter[1] = 1
    Counter[2] = 2
    Counter[3] = 3
    Counter[4] = 4

    general_library = {}
    general_library['timestamp'] = timestamp
    general_library['ME_101'] = ME_101
    general_library['TE_101'] = TE_101
    general_library['TE_102'] = TE_102
    general_library['PT_101'] = PT_101
    general_library['PT_102'] = PT_102
    general_library['PT_103'] = PT_103
    general_library['MFC_101'] = MFC_101
    general_library['MFC_102'] = MFC_102
    general_library['MFC_103'] = MFC_103
    general_library['MFC_104'] = MFC_104
    general_library['FE_101'] = FE_101
    general_library['CO_NDIR'] = CO_NDIR
    general_library['CO2_NDIR'] = CO2_NDIR
    general_library['CH4_NDIR'] = CH4_NDIR
    general_library['H2_GC'] = H2_GC
    general_library['CO_GC'] = CO_GC
    general_library['CO2_GC'] = CO2_GC
    general_library['CH4_GC'] = CH4_GC
    general_library['C2H6_GC'] = C2H6_GC
    general_library['N2_GC'] = N2_GC
    general_library['Ar_GC'] = Ar_GC
    general_library['counter'] = Counter

    units = {}
    units['timestamp'] = None
    units['ME_101'] = 'lb/hr'
    units['TE_101'] = 'K'
    units['TE_102'] = 'K'
    units['PT_101'] = 'lbf/in^2'
    units['PT_102'] = 'lbf/in^2'
    units['PT_103'] = 'lbf/in^2'
    units['MFC_101'] = 'L/min'
    units['MFC_102'] = 'L/min'
    units['MFC_103'] = 'L/min'
    units['MFC_104'] = 'L/min'
    units['FE_101'] = 'L/min'
    units['CO_NDIR'] = None
    units['CO2_NDIR'] = None
    units['CH4_NDIR'] = None
    units['H2_GC'] = None
    units['CO_GC'] = None
    units['CO2_GC'] = None
    units['CH4_GC'] = None
    units['C2H6_GC'] = None
    units['N2_GC'] = None
    units['Ar_GC'] = None
    units['counter'] = None

    start = datetime.datetime(1981,07,06,13,13,12)
    end = datetime.datetime(1981,07,06,13,13,16)
    
    """Needs to:
    + Load a time series data frame correctly
        - Hold the correct data
        - Hold the correct data types
        - Hold no extraneous data
        - Ensure that the counter is appropriately initialized
    + Raise an exception on an improper interface (no type)
    + Raise an exception on an unconnected interface
    + Raise an exception on an interface that has no database selected
    + Raise an exception on a bad table name
    + Raise an exception if there is not a "timestamp" field in the database
    + Raise an exception if the start time makes no sense
    + Raise an exception if the stop time makes no sense
    + Raise an exception if start time > stop time
    
    """
    def testDataLoadCorrect(self):
        """Test whether all the data is loaded correctly and no extraneous data exists"""
        interface = db.db_interface(host = "192.168.10.20", user = "ryon", passwd = "3lectron*")
        interface.connect()
        q = db.use_Query("gas_unit_test")
        interface.query(q)
        
        ts = lfl.ts_data(start = LoadDataTests.start, end = LoadDataTests.end)
        ts.SQL_load(interface, table = "gas_data_test_tbl")
        for key in LoadDataTests.general_library.keys():
                for i in range(len(ts[key])):
                    self.assertEqual(ts[key][i], LoadDataTests.general_library[key][i])
        for key in ts.data.keys():
            self.assertIn(key,LoadDataTests.general_library.keys())

    def testBadInterface(self):
        """Test whether an interface is bad or if it is connected"""
        interface = ""
        ts = lfl.ts_data(start = LoadDataTests.start, end = LoadDataTests.end)
        self.assertRaises(lfl.SQLInterfaceError, ts.SQL_load,interface, table = "gasifier_lv_GC_view")

    def testUnconnectedInterface(self):
        """An unconnected interface should raise an error"""
        interface = db.db_interface(host = "192.168.13.15", user = "chris", passwd = "udflyer87")
        ts = lfl.ts_data(start = LoadDataTests.start, end = LoadDataTests.end)
        self.assertRaises(lfl.SQLInterfaceError, ts.SQL_load,interface, table = "gasifier_lv_GC_view")

    def testNoDatabaseSelected(self):
        """An interface with no selected database should raise an error:"""
        interface = db.db_interface(host = "192.168.13.15", user = "chris", passwd = "udflyer87")
        interface.connect()
        ts = lfl.ts_data(start = LoadDataTests.start, end = LoadDataTests.end)
        self.assertRaises(lfl.SQLInterfaceError, ts.SQL_load,interface, table = "gasifier_lv_GC_view")

    def testDataIsTimeseries(self):
        """The data should have a timestamp column (actually be timeseries data)"""
        interface = db.db_interface(host = "192.168.10.20", user = "ryon", passwd = "3lectron*")
        interface.connect()
        q = db.use_Query("gas_unit_test")
        interface.query(q)
        ts = lfl.ts_data(start = LoadDataTests.start, end = LoadDataTests.end)
        ts.SQL_load(interface, table = "gas_data_test_tbl")
        self.assertIn("timestamp", ts.data.keys())

    def testSensibleStartandEnd(self):
        """The start and end times should be datetime.datetime objects"""
        
        self.assertRaises(lfl.TimeError, lfl.ts_data, "fart", LoadDataTests.end)
        self.assertRaises(lfl.TimeError, lfl.ts_data, LoadDataTests.start, end = "fart")

    def testChronologicaltimes(self):
        """The start should happen before the end"""
        
        self.assertRaises(lfl.TimeError, lfl.ts_data, LoadDataTests.end, LoadDataTests.start)
        
   

class InterpolateDataTests(unittest.TestCase):
    y = np.array([2.1,3.2,None,None,None,6.5,7.2,None,None,None,3.5],dtype = 'float64')
    x = np.array([1,2,3,4,5,6,7,8,9,11,12],dtype = 'float64')
    leadingNone = np.array([None,3.2,None,None,None,6.5,7.2,None,None,None,3.5],dtype = 'float64')
    trailingNone = np.array([2.1,3.2,None,None,None,6.5,7.2,None,None,None,None],dtype = 'float64')
    answer_y = np.array([2.1,3.2,4.025,4.85,5.675,6.5,7.2,6.46,5.72,4.24,3.5],dtype = 'float64')
    answer_leadN = np.array([None,3.2,4.025,4.85,5.675,6.5,7.2,6.46,5.72,4.24,3.5],dtype = 'float64')
    answer_trailN = np.array([2.1,3.2,4.025,4.85,5.675,6.5,7.2,None,None,None,None],dtype = 'float64')
    start = datetime.datetime(1981,07,06,13,12,12)
    end = datetime.datetime(1981,07,06,13,12,16)

    """Needs to:
    + Correctly interpolate within a series column with interspersed None
        + Give the correct interpolated values at different positions
        + Keep leading and trailing None values
    + Raise an exception if the argument is not a data column contained in the time series
    + Do nothing to the data if there are no None's in it
    """

    def testCorrectInterpolationNoLeadorTrail(self):
        """The returned interpolated column must be correctly calculated"""
        ts = lfl.ts_data(start = InterpolateDataTests.start, end = InterpolateDataTests.end)
        ts['x'] = InterpolateDataTests.x
        ts['y'] = InterpolateDataTests.y
        for i in range(len(InterpolateDataTests.answer_y)):

            self.assertAlmostEqual(ts.interpolate_col(x = "x", y = "y")[i], InterpolateDataTests.answer_y[i], 4)
        


    def testCorrectInterpolationLead(self):
        """The returned interpolated column must be correctly calculated with a leading set of Nones"""
        ts = lfl.ts_data(start = InterpolateDataTests.start, end = InterpolateDataTests.end)
        ts['x'] = InterpolateDataTests.x
        ts['y'] = InterpolateDataTests.leadingNone
        for i in range(len(InterpolateDataTests.answer_leadN)):
           if np.isnan(InterpolateDataTests.answer_leadN[i]):
               self.assertTrue(np.isnan(ts.interpolate_col(x="x",y="y")[i]))
           else:
               self.assertAlmostEqual(ts.interpolate_col(x = "x", y = "y")[i], InterpolateDataTests.answer_leadN[i], 4)

    def testCorrectInterpolationTrail(self):
        """The returned interpolated column must be correctly calculated with a trailing set of Nones"""
        ts = lfl.ts_data(start = InterpolateDataTests.start, end = InterpolateDataTests.end)
        ts['x'] = InterpolateDataTests.x
        ts['y'] = InterpolateDataTests.trailingNone
        for i in range(len(InterpolateDataTests.answer_trailN)):
            if np.isnan(InterpolateDataTests.answer_trailN[i]):
                self.assertTrue(np.isnan(ts.interpolate_col(x="x",y="y")[i]))
            else:
                self.assertEqual(ts.interpolate_col(x = "x", y = "y")[i], InterpolateDataTests.answer_trailN[i])

    def testCorrectInterpolationNothing(self):
        """The returned interpolated column must be correctly calculated when there is nothing to do"""
        ts = lfl.ts_data(start = InterpolateDataTests.start, end = InterpolateDataTests.end)
        ts['x'] = InterpolateDataTests.x
        ts['y'] = InterpolateDataTests.answer_y
        for i in range(len(InterpolateDataTests.answer_y)):
            self.assertEqual(ts.interpolate_col(x = "x", y = "y")[i], InterpolateDataTests.answer_y[i])

    def testBadArgument(self):
        """An exception must be raised on an invalid column name"""
        ts = lfl.ts_data(start = InterpolateDataTests.start, end = InterpolateDataTests.end)
        ts['x'] = InterpolateDataTests.x
        ts['y'] = InterpolateDataTests.answer_y
        self.assertRaises(lfl.InterpArgumentError,ts.interpolate_col,x = "x", y = "noggin")
    
#This is NOT for time series data -- this is a property of the more general analyzer
class LoadRunParametersTests(unittest.TestCase):
    """Needs to:
    0 Correctly load a set of run parameters into a dictionary from the database
        - Hold the correct parameters
        - Hold nothing else than the correct parameters
    0 Raise an exception on an improper interface (no type)
    0 Raise an exception on an unconnected interface
    0 Raise an exception on an interface that has no database selected
    0 Raise an exception on a bad table name
    """
    pass

class GlossaryReplaceTests(unittest.TestCase):
    """Needs to:
    + Correctly replace the names of tags with glossary names
        - Don't change any names that are NOT in the glossary
        - Change all names that ARE in the glossary
    + Raise an exception if a tag in the glossary does not match a tag in the time-series
    + Raise an exception if the glossary is not a dict 
    """
    
    glossary = {"TE_101": "inlet_temp", "TE_102":"outlet_temp", "PT_102":"outlet_pressure"}
    pre_keys = ("TE_101", "TE_102", "PT_102")
    post_keys = ("inlet_temp", "outlet_temp","PT_101", "outlet_pressure")

    def testCorrectGlossaryReplace(self):
        """The names (and only the names) in the glossary should be changed"""
        ts = lfl.ts_data(start=datetime.datetime(1981,07,06,13,12,12), end = datetime.datetime(1981,07,06,13,12,16), array_dict = LoadDataTests.general_library.copy())
        ts.glossary_replace(GlossaryReplaceTests.glossary)
        for key in GlossaryReplaceTests.post_keys:
            self.assertIn(key,ts.data.keys())
        for key in GlossaryReplaceTests.pre_keys:
            self.assertNotIn(key, ts.data.keys())

    def testBadGlossaryName(self):
        """A tag in the glossary NOT in the timeseries should raise an error"""
        ts = lfl.ts_data(start=datetime.datetime(1981,07,06,13,12,12), end = datetime.datetime(1981,07,06,13,12,16), array_dict = LoadDataTests.general_library.copy())
        glossary1 = {"TE_101":"inlet_temp", "HI_203":"help_indicator"}
        self.assertRaises(df.GlossaryError, ts.glossary_replace, glossary1)

    def testBadGlossaryType(self):
        """If the glossary is not a dictionary, raise an error"""
        ts = lfl.ts_data(start=datetime.datetime(1981,07,06,13,12,12), end = datetime.datetime(1981,07,06,13,12,16), array_dict = LoadDataTests.general_library.copy())
        glossary = 2.4
        self.assertRaises(df.GlossaryError, ts.glossary_replace, glossary)
        

class FeedrateUnitConversionTests(unittest.TestCase):
    """Needs to:
    + Correctly convert units on a set of data in a defined timeseries
        - Make factually correct conversions and return the correct column
        - Make sure that the units are now changed to the new unit
        
    + Raise an exception if the time series does not have a given column
    + Raise an exception on inlet conversions that are not in the conversion dictionary
    + Raise an exception on inconsistent units
    """
    def testCorrectUnitConversion(self):
        """Unit conversion must be performed correctly, column must be added to the series, and the new units must be correct"""
        ts = lfl.ts_data(start=datetime.datetime(1981,07,06,13,12,12), end = datetime.datetime(1981,07,06,13,12,16),array_dict = LoadDataTests.general_library, units_dict = LoadDataTests.units.copy())
        for key, value in LoadDataTests.general_library.items():
            ts[key] = value
                
        lbh_kgs_conversion = LoadDataTests.general_library['ME_101'].copy()/2.20462/3600.0  #this is the correct answer
        ts.convert_col_units('ME_101', 'kg/s')
        for i in range(len(lbh_kgs_conversion)):
            self.assertEqual(ts['ME_101'][i],lbh_kgs_conversion[i])
        self.assertEqual(ts.units['ME_101'], 'kg/s')
    
    def testNoColumnError(self):
        """Trying to convert a non-existent column should raise an error"""
        ts = lfl.ts_data(start=datetime.datetime(1981,07,06,13,12,12), end = datetime.datetime(1981,07,06,13,12,16),array_dict = LoadDataTests.general_library, units_dict = LoadDataTests.units.copy())
        for key, value in LoadDataTests.general_library.items():
            ts[key] = value
        
        
        self.assertRaises(df.NoColumnError, ts.convert_col_units, 'ME_520','kg/s')

    def testNoConversionInformation(self):
        """Trying to convert units that do not have a conversion factor should raise an error"""
        ts = lfl.ts_data(start=datetime.datetime(1981,07,06,13,12,12), end = datetime.datetime(1981,07,06,13,12,16),array_dict = LoadDataTests.general_library, units_dict = LoadDataTests.units.copy())
        for key, value in LoadDataTests.general_library.items():
            ts[key] = value
        
        self.assertRaises(lfl.UnitConversionError, ts.convert_col_units, 'ME_101', 'stone/s')


    def testInconsistentUnitConversion(self):
        """Trying to convert to inconsistent units should raise an error"""
        ts = lfl.ts_data(start=datetime.datetime(1981,07,06,13,12,12), end = datetime.datetime(1981,07,06,13,12,16),array_dict = LoadDataTests.general_library, units_dict = LoadDataTests.units.copy())

        for key, value in LoadDataTests.general_library.items():
            ts[key] = value
        
        self.assertRaises(lfl.UnitConversionError, ts.convert_col_units, 'ME_101', 'm/s')

class ElementalFeedrateCalculationTests(unittest.TestCase):
    """Needs to:
    + Correctly calculate the feed rate of the given element in the inlet stream
    + Raise an exception on an element not in the parameter list
    + Raise an exception on a bad argument list
    0 Raise an exception on an empty stream list
    """
    def testProperElementalFlowrates(self):
        """Test that the elemental flowrates in and out are properly calculated"""
        pts = lfl.ProcTS(start = datetime.datetime(1981,07,06,13,12,12), end = datetime.datetime(1981,07,06,13,12,16))
        for key,value in LoadDataTests.general_library.items():
            pts[key] = value
        pts.set_units(LoadDataTests.units)
        biomass_feed = lfl.Stream('ME_101',flowrate = pts.get_val('ME_101'), composition = {'biomass':1.00}, basis = "mass")
        gas_exit = lfl.Stream('FE_101', flowrate = pts.get_val('FE_101'), basis = "gas_volume")
        biomass_feed.temperature = (25, "C")
        gas_exit.temperature = (25, "C")
        biomass_feed.pressure = (101325, 'kg/s^2/m')
        gas_exit.pressure = (101325, 'kg/s^2/m')

        exit_list = ['H2', 'CO', 'CO2', 'CH4', 'C2H6', 'N2', 'Ar']
        for specie in exit_list:
            key = "%s_GC" % specie
            gas_exit.composition[specie] = pts[key]
        pts.inlet_streams = [biomass_feed]
        pts.outlet_streams = [gas_exit]

        biomass_breakdown = {}
        biomass_breakdown['biomass'] = {'C':0.5, 'O': 0.4, 'H': 0.1}
        biomass_feed.special_species = biomass_breakdown
        pts.proc_elements = ['C', 'H', 'O']
        pts.generate_inlet_outlet_elemental_flows()
        biomass_molar_C = np.array([0.018379,0.0178546,0.018905,0.01838,0.018379875]) #Need to define this answer
        gas_molar_C = np.array([0.01431, 0.01468, 0.014287, 0.014536, 0.0151035])
        gas_molar_O = np.array([0.015288, 0.015657, 0.015119, 0.0155047, 0.0161151])
        gas_molar_H = np.array([0.036413, 0.035788, 0.037243, 0.036824, 0.037183])
        for i in range(len(biomass_molar_C)):
            self.assertAlmostEqual(pts['C_inlet'][i],biomass_molar_C[i], places = 4)
            self.assertAlmostEqual(pts['C_outlet'][i],gas_molar_C[i], places = 4)
            self.assertAlmostEqual(pts['O_outlet'][i],gas_molar_O[i], places = 4)
            self.assertAlmostEqual(pts['H_outlet'][i],gas_molar_H[i], places = 4)

    def testBadStreamInput(self):
        """An improper stream should throw an error"""
        pts = lfl.ProcTS(start = datetime.datetime(1981,07,06,13,12,12), end = datetime.datetime(1981,07,06,13,12,16))
        for key,value in LoadDataTests.general_library.items():
            pts[key] = value
        pts.set_units(LoadDataTests.units)
        biomass_feed = lfl.Stream('ME_101', flowrate = pts.get_val('ME_101'), composition = {'biomass':1.00}, basis = "mass")
        gas_exit = lfl.Stream('FE_101', flowrate = pts.get_val('FE_101'), basis = "gas_volume")
        exit_list = ['H2', 'CO', 'CO2', 'CH4', 'C2H6', 'N2', 'Ar']
        for specie in exit_list:
            key = "%s_GC" % specie
            gas_exit.composition[specie] = pts[key]
        pts.proc_elements = ['C', 'H', 'O']
        pts.inlet_streams = ["a"]
        pts.outlet_streams = [gas_exit]
        biomass_breakdown = {}
        biomass_breakdown['biomass'] = {'C':0.5, 'O': 0.4, 'H': 0.1}
        biomass_feed.special_species = biomass_breakdown
        self.assertRaises(lfl.BadStreamError, pts.generate_inlet_outlet_elemental_flows)
       
        pts.inlet_streams = [biomass_feed]
        pts.outlet_streams = ["a"]
        self.assertRaises(lfl.BadStreamError, pts.generate_inlet_outlet_elemental_flows)

    def testBadElementInput(self):
        """A non-physical element should throw an error"""
        pts = lfl.ProcTS(start = datetime.datetime(1981,07,06,13,12,12), end = datetime.datetime(1981,07,06,13,12,16))
        for key,value in LoadDataTests.general_library.items():
            pts[key] = value
        pts.set_units(LoadDataTests.units)
        biomass_feed = lfl.Stream('ME_101',flowrate = pts.get_val('ME_101'), composition = {'biomass':1.00}, basis = "mass")
        gas_exit = lfl.Stream('FE_101', flowrate = pts.get_val('FE_101'), basis = "gas_volume")
        exit_list = ['H2', 'CO', 'CO2', 'CH4', 'C2H6', 'N2', 'Ar']
        for specie in exit_list:
            key = "%s_GC" % specie
            gas_exit.composition[specie] = pts[key]
        pts.inlet_streams = [biomass_feed]
        pts.outlet_streams = [gas_exit]
        pts.proc_elements = ['Z', 'O', 'H']
        self.assertRaises(lfl.BadStreamError, pts.generate_inlet_outlet_elemental_flows)

    def testStreamsEmpty(self):
        """An empty stream list on the inlet or outlet should throw an error"""
        pass

class InletMolarFlowrateTests(unittest.TestCase):
    """Needs to:
    + Correctly calculate the molar flowrate of all species in a given time series list
        + Use argument dictionary to determine which tags contain which components
        + Correctly calculate the components from all streams as molar flowrates (mol/s)
        - Correctly add the proper units to the time series dataframe
    
    + Raise exception on incorrect argument types
    
    """
    def testProperMolarFlowrates(self):
        
        """Test that the species flowrates in and out are properly calculated"""
        pts = lfl.ProcTS(start = datetime.datetime(1981,07,06,13,12,12), end = datetime.datetime(1981,07,06,13,12,16))
        for key,value in LoadDataTests.general_library.items():
            pts[key] = value
        pts.set_units(LoadDataTests.units)
        entrainment_gas_feed = lfl.Stream('MFC_101', flowrate = pts.get_val('MFC_101'), composition = {'N2':0.95, 'Ar':0.05}, basis = "gas_volume")
        methane_gas_feed = lfl.Stream('MFC_102', flowrate = pts.get_val('MFC_102'), composition = {'CH4':1.00}, basis = "gas_volume")
        gas_exit = lfl.Stream('FE_101', flowrate = pts.get_val('FE_101'), basis = "gas_volume")
        entrainment_gas_feed.temperature = (25,"C")
        methane_gas_feed.temperature = (25, "C")
        gas_exit.temperature = (25, "C")
        entrainment_gas_feed.pressure = (101325, "Pa")
        methane_gas_feed.pressure = (101325, "Pa")
        gas_exit.pressure = (101325, 'kg/s^2/m')
        
        exit_list = ['H2', 'CO', 'CO2', 'CH4', 'C2H6', 'N2', 'Ar']
        for specie in exit_list:
            key = "%s_GC" % specie
            gas_exit.composition[specie] = pts[key]
        pts.inlet_streams = [entrainment_gas_feed, methane_gas_feed]
        pts.outlet_streams = [gas_exit]

        
        pts.proc_elements = ['C', 'H', 'O']
        pts.proc_species = ['H2', 'CO', 'CO2', 'CH4', 'C2H6', 'N2']

        pts.generate_inlet_outlet_species_flows()
        inlet_species_CH4 = np.array([0.0048370,0.004830218,0.00484384,0.0048370,0.0048370]) 
        inlet_species_CO2 = np.array([0,0,0,0,0]) 
        inlet_species_CO = np.array([0,0,0,0,0]) 
        inlet_species_N2 = np.array([0,0,0,0,0]) 
        outlet_species_N2 = LoadDataTests.general_library['FE_101']*LoadDataTests.general_library['N2_GC']/60000.0*101325/(273.15+25.0)/8.314 
        outlet_species_CO = LoadDataTests.general_library['FE_101']*LoadDataTests.general_library['CO_GC']/60000.0*101325/(273.15+25.0)/8.314 
        outlet_species_CO2 = LoadDataTests.general_library['FE_101']*LoadDataTests.general_library['CO2_GC']/60000.0*101325/(273.15+25.0)/8.314 
        outlet_species_CH4 = LoadDataTests.general_library['FE_101']*LoadDataTests.general_library['CH4_GC']/60000.0*101325/(273.15+25.0)/8.314
        for i in range(len(inlet_species_CH4)):
 
            self.assertAlmostEqual(pts['CH4_inlet'][i],inlet_species_CH4[i])
            self.assertAlmostEqual(pts['CO2_inlet'][i],inlet_species_CO2[i])
            self.assertAlmostEqual(pts['CO_inlet'][i],inlet_species_CO[i])
            self.assertAlmostEqual(pts['N2_inlet'][i],inlet_species_N2[i])
            self.assertAlmostEqual(pts['CH4_outlet'][i],outlet_species_CH4[i])
            self.assertAlmostEqual(pts['CO2_outlet'][i],outlet_species_CO2[i])
            self.assertAlmostEqual(pts['CO_outlet'][i],outlet_species_CO[i])
            self.assertAlmostEqual(pts['N2_outlet'][i],outlet_species_N2[i])

   

    def testBadStream(self):
        pts = lfl.ProcTS(start = datetime.datetime(1981,07,06,13,12,12), end = datetime.datetime(1981,07,06,13,12,16))
        for key,value in LoadDataTests.general_library.items():
            pts[key] = value
        pts.set_units(LoadDataTests.units)
        entrainment_gas_feed = lfl.Stream('MFC_101', flowrate = pts.get_val('MFC_101'),composition = {'N2':0.95, 'Ar':0.05}, basis = "gas_volume")
        methane_gas_feed = lfl.Stream('MFC_102', flowrate = pts.get_val('MFC_102'),composition = {'CH4':1.00}, basis = "gas_volume")
        gas_exit = lfl.Stream('FE_101', flowrate = pts.get_val('FE_101'),basis = "gas_volume")
        entrainment_gas_feed.temperature = (25,"C")
        methane_gas_feed.temperature = (25, "C")
        gas_exit.temperature = (25, "C")
        entrainment_gas_feed.pressure = (101325, "Pa")
        methane_gas_feed.pressure = (101325, "Pa")
        gas_exit.pressure = (101325, 'kg/s^2/m')
        
        exit_list = ['H2', 'CO', 'CO2', 'CH4', 'C2H6', 'N2', 'Ar']
        for specie in exit_list:
            key = "%s_GC" % specie
            gas_exit.composition[specie] = pts[key]
               
        pts.proc_elements = ['C', 'H', 'O']
        pts.proc_species = ['H2', 'CO', 'CO2', 'CH4', 'C2H6', 'N2']

        
        
        pts.inlet_streams = ['fart', methane_gas_feed]
        pts.outlet_streams = [gas_exit]
        self.assertRaises(lfl.BadStreamError, pts.generate_inlet_outlet_species_flows)

        pts.inlet_streams = [entrainment_gas_feed, methane_gas_feed]
        pts.outlet_streams = ['hello']
        self.assertRaises(lfl.BadStreamError, pts.generate_inlet_outlet_species_flows)

    def TestProperUnits(self):
        """All the units in the process stream must be mol/s"""
        pts = lfl.ProcTS(start = datetime.datetime(1981,07,06,13,12,12), end = datetime.datetime(1981,07,06,13,12,16))
        for key,value in LoadDataTests.general_library.items():
            pts[key] = value
        pts.set_units(LoadDataTests.units)
        entrainment_gas_feed = lfl.Stream('MFC_101', flowrate = pts.get_val('MFC_101'),composition = {'N2':0.95, 'Ar':0.05}, basis = "gas_volume")
        methane_gas_feed = lfl.Stream('MFC_102', flowrate = pts.get_val('MFC_102'),composition = {'CH4':1.00}, basis = "gas_volume")
        gas_exit = lfl.Stream('FE_101', flowrate = pts.get_val('FE_101'), basis = "gas_volume")
        entrainment_gas_feed.temperature = (25,"C")
        methane_gas_feed.temperature = (25, "C")
        gas_exit.temperature = (25, "C")
        entrainment_gas_feed.pressure = (101325, "Pa")
        methane_gas_feed.pressure = (101325, "Pa")
        gas_exit.pressure = (101325, 'kg/s^2/m')
        
        exit_list = ['H2', 'CO', 'CO2', 'CH4', 'C2H6', 'N2', 'Ar']
        for specie in exit_list:
            key = "%s_GC" % specie
            gas_exit.composition[specie] = pts[key]
        pts.inlet_streams = [entrainment_gas_feed, methane_gas_feed]
        pts.outlet_streams = [gas_exit]

        
        pts.proc_elements = ['C', 'H', 'O']
        pts.proc_species = ['H2', 'CO', 'CO2', 'CH4', 'C2H6', 'N2']

        pts.generate_inlet_outlet_species_flows()
        
        self.assertEqual(pts.units['CH4_inlet'],'mol/s')
        self.assertEqual(pts.units['CO2_inlet'], 'mol/s')
        self.assertEqual(pts.units['CO_inlet'], 'mol/s')
        self.assertEqual(pts.units['N2_inlet'], 'mol/s')
        self.assertEqual(pts.units['CH4_outlet'], 'mol/s')
        self.assertEqual(pts.units['CO2_outlet'], 'mol/s')
        self.assertEqual(pts.units['CO_outlet'], 'mol/s')
        self.assertEqual(pts.units['N2_outlet'], 'mol/s')    

class ConversionTests(unittest.TestCase):
    """Needs to:
    + Correctly calculate all types of conversion given fixed inputs:
        - Total conversion
        - Good conversion
        - CO yield
        - Standard conversion
    0 Raise an exception on bad argument types
    
    + Raise an exception if the conversion function is called before the feedrates are set (elemental C, species)
    """
    def testCorrectCarbonConversions(self):
        """Carbon conversions must be calculated accurately and appropriate units set"""
        
        gts = lfl.GasifierProcTS(start = datetime.datetime(1981,07,06,13,12,12),end=datetime.datetime(1981,07,06,13,12,16))
        for (key, value) in LoadDataTests.general_library.items():
            gts[key] = value
        gts.set_units(LoadDataTests.units)
        biomass_feed = lfl.Stream('ME_101',flowrate = gts.get_val('ME_101'), composition = {'biomass':1.00}, basis = "mass")
        entrainment_gas_feed = lfl.Stream('MFC_101', flowrate = gts.get_val('MFC_101'), composition = {'N2':0.95, 'Ar':0.05}, basis = "gas_volume")
        methane_gas_feed = lfl.Stream('MFC_102', flowrate = gts.get_val('MFC_102'), composition = {'CH4':1.00}, basis = "gas_volume")
        gas_exit = lfl.Stream('FE_101', flowrate = gts.get_val('FE_101'),basis = "gas_volume")
        
        biomass_feed.temperature = (25,"C")
        biomass_feed.pressure = (101325, "Pa")


        entrainment_gas_feed.temperature = (25,"C")
        methane_gas_feed.temperature = (25, "C")
        gas_exit.temperature = (25, "C")
        entrainment_gas_feed.pressure = (101325, "Pa")
        methane_gas_feed.pressure = (101325, "Pa")
        gas_exit.pressure = (101325, 'kg/s^2/m')

        exit_list = ['H2', 'CO', 'CO2', 'CH4', 'C2H6', 'N2', 'Ar']
        for specie in exit_list:
            key = "%s_GC" % specie
            gas_exit.composition[specie] = gts[key]
        gts.inlet_streams = [entrainment_gas_feed, methane_gas_feed, biomass_feed]
        gts.outlet_streams = [gas_exit]
        gts.proc_elements = ['C', 'H', 'O']
        gts.proc_species = ['H2', 'CO', 'CO2', 'CH4', 'C2H6', 'N2']
        biomass_breakdown = {}
        biomass_breakdown['biomass'] = {'C':0.5, 'O': 0.4, 'H': 0.1}
        biomass_feed.special_species = biomass_breakdown
        gts.generate_inlet_outlet_elemental_flows()
        gts.generate_inlet_outlet_species_flows()
        
     
        gts.generate_carbon_conversions()
        inlet = LoadDataTests.general_library['ME_101']*0.5*453.592909/3600.0/12.0 + LoadDataTests.general_library['MFC_102']/60000.0*101325/8.314/298.15
        outlet_molar = LoadDataTests.general_library['FE_101']/60000*101325/8.314/298.15
        
        goodX = outlet_molar*(LoadDataTests.general_library['CO_GC']+LoadDataTests.general_library['CO2_GC'])/inlet 
        totX = outlet_molar*(LoadDataTests.general_library['CO_GC']+LoadDataTests.general_library['CO2_GC']+LoadDataTests.general_library['CH4_GC']+2.0*LoadDataTests.general_library['C2H6_GC'])/inlet 
        stdX = outlet_molar*(LoadDataTests.general_library['CO_GC']+LoadDataTests.general_library['CO2_GC']+LoadDataTests.general_library['CH4_GC'])/inlet 
        phiCO = outlet_molar*(LoadDataTests.general_library['CO_GC'])/inlet

        for i in range(len(goodX)):
            self.assertAlmostEqual(gts['X_good'][i],goodX[i],3)
            self.assertAlmostEqual(gts['X_tot'][i],totX[i],3)
            self.assertAlmostEqual(gts['X_std'][i],stdX[i],3)
            self.assertAlmostEqual(gts['CO_yield'][i],phiCO[i],3)

        for name in ['X_good','X_tot','X_std','CO_yield']:
            self.assertEqual(gts.units[name],None)

    def testConversionBeforeInletsCalc(self):
        """Trying to calculate conversions before generating the proper flows should throw an error"""
        gts = lfl.GasifierProcTS(start = datetime.datetime(1981,07,06,13,12,12),end=datetime.datetime(1981,07,06,13,12,16))
        for (key, value) in LoadDataTests.general_library.items():
            gts[key] = value
        gts.set_units(LoadDataTests.units)
        biomass_feed = lfl.Stream('ME_101',flowrate = gts.get_val('ME_101'),composition = {'biomass':1.00}, basis = "mass")
        entrainment_gas_feed = lfl.Stream('MFC_101', flowrate = gts.get_val('MFC_101'),composition = {'N2':0.95, 'Ar':0.05}, basis = "gas_volume")
        methane_gas_feed = lfl.Stream('MFC_102', flowrate = gts.get_val('MFC_102'),composition = {'CH4':1.00}, basis = "gas_volume")
        gas_exit = lfl.Stream('FE_101', flowrate = gts.get_val('FE_101'),basis = "gas_volume")
        
        biomass_feed.temperature = (25,"C")
        biomass_feed.pressure = (101325, "Pa")


        entrainment_gas_feed.temperature = (25,"C")
        methane_gas_feed.temperature = (25, "C")
        gas_exit.temperature = (25, "C")
        entrainment_gas_feed.pressure = (101325, "Pa")
        methane_gas_feed.pressure = (101325, "Pa")
        gas_exit.pressure = (101325, 'kg/s^2/m')

        exit_list = ['H2', 'CO', 'CO2', 'CH4', 'C2H6', 'N2', 'Ar']
        for specie in exit_list:
            key = "%s_GC" % specie
            gas_exit.composition[specie] = gts[key]
        gts.inlet_streams = [entrainment_gas_feed, methane_gas_feed, biomass_feed]
        gts.outlet_streams = [gas_exit]
        gts.proc_elements = ['C', 'H', 'O']
        gts.proc_species = ['H2', 'CO', 'CO2', 'CH4', 'C2H6', 'N2']
        biomass_breakdown = {}
        biomass_breakdown['biomass'] = {'C':0.5, 'O': 0.4, 'H': 0.1}
        biomass_feed.special_species = biomass_breakdown
        
        self.assertRaises(lfl.ConversionError, gts.generate_carbon_conversions)
        gts.generate_inlet_outlet_elemental_flows()
        self.assertRaises(lfl.ConversionError, gts.generate_carbon_conversions)
     
        
        

class DataFilteringTests(unittest.TestCase):
    """Needs to:
    + Correctly filter arbitrary low value off a column to Nan
    + Correctly filter arbitrary high value off a column to Nan
    + Raise an exception on an invalid argument type
    + Raise an exception on a non-existent column
    """

    def testLowFilter(self):
        """Must properly filter a low value from a data frame"""
        a = np.array([2.5,3.6,100.2,3.5])
        b = np.array([-3.2,0.0, -2.3, 1.5])
        data1 = df.Dataframe(array_dict = {'a':a, 'b':b})
        data1['b'] = data1.filter_vals(column = 'b', v = 0, action = 'low')
        for i in [1,3]:
            self.assertEqual(data1['b'][i], np.array([np.nan, 0, np.nan, 1.5])[i])

        self.assertTrue(np.isnan(data1['b'][0]))
	self.assertTrue(np.isnan(data1['b'][2]))

    def testHighFilter(self):
        """Must properly filter a high value from a data frame"""
        a = np.array([2.5,3.6,100.2,3.5])
        b = np.array([-3.2,0.0, -2.3, 1.5])
        data1 = df.Dataframe(array_dict = {'a':a, 'b':b})
        data1['a'] = data1.filter_vals(column = 'a', v = 3.5, action = 'high')
        for i in [0,3]:
            self.assertEqual(data1['a'][i], np.array([2.5, np.nan, np.nan, 3.5])[i])
        self.assertTrue(np.isnan(data1['a'][1]))
        self.assertTrue(np.isnan(data1['a'][2]))

    def testBadArgumentType(self):
        """A non-numeric filter value should raise an error or non-string field name"""
        a = np.array([2.5,3.6,100.2,3.5])
        b = np.array([-3.2,0.0, -2.3, 1.5])
        data1 = df.Dataframe(array_dict = {'a':a, 'b':b})
        
        self.assertRaises(df.BadArgumentError, data1.filter_vals, 'a', 'hi', 'high')
        self.assertRaises(df.BadArgumentError, data1.filter_vals, 2, 3.5, 'high')
        self.assertRaises(df.BadArgumentError, data1.filter_vals, 'a', 3.5, 23.5)
        self.assertRaises(df.BadArgumentError, data1.filter_vals, 'a',3.5,'mid')
    
    def testNoColumn(self):
        """A column that is not in the data frame should raise an error"""

        a = np.array([2.5,3.6,100.2,3.5])
        b = np.array([-3.2,0.0, -2.3, 1.5])
        data1 = df.Dataframe(array_dict = {'a':a, 'b':b})
        
        self.assertRaises(df.NoColumnError, data1.filter_vals, 'c', 2.5, 'high')

    

class CompositionTests(unittest.TestCase):
    """Needs to:
    0 Correctly normalize a given column with an identified set of inert flowrates
    0 Raise an exception for a column key not in the time series
    0 Raise an exception on divide by zero
    """
    def testCorrectComposition(self):
        """Normalized compositions must be correctly calculated"""
        pass

    def testBadColumn(self):
        """If a "excluded" column or a "normalized" column is not in the dataframe, raise an error"""
        pass

    def testDivByZero(self):
        """If the denominator in the normalization column is zero, raise an error"""
        pass


class EnthalpyTests(unittest.TestCase):
    """Needs to:
    0 Correctly calculate sensible heat change of a set of marked species given:
        - Inlet temperature
        - Outlet temperature
        - Vector of desired species
        - Object containing Cp for all species (dictionary of NASA polynomials?)
    0 Raise exception if outlet/inlet molar gas flowrates not yet calculated
    0 Raise exception on non-sensical temperature arguments
    0 Raise exception if species flowrate not there
    0 Raise exception on non-sensical species vector
    0 Raise exception on wrong thermo data type (can we directly integrate Cantera??? - if we can, don't need to do this)
    0 Correctly calculate pressure heat change of a set of marked species given:
        - Inlet pressure
        - Constant temperature
        - Outlet pressure
        - Vector of desired species
    0 Raise exception on non-sensical pressure arguments
    0 Raise exception if species flowrate not there
    0 Raise exception on non-sensical species vector
    0 Correctly calculate overall enthalpy change from reactions given:
        - Inlet temperature
        - Outlet temperature
        - Inlet Pressure
        - Outlet Pressure
        - Vector of desired species
        - Object containing Cp, dH_formation for all species
    0 Raise exception on non-sensical arguments
    """

    def testStreamEnthalpy(self):
        """Must correctly calculate stream enthalpy"""
        test_stream = lfl.Stream('test_stream', flowrate = (10, 'L/min'), temperature = (400, 'K'), pressure = (101325, 'Pa'), composition = {'CO2':0.5, 'H2O':0.5}, basis = 'std_gas_volume', std_temperature = (70, 'F'))
        
        #Hand calculate enthlapy (kJ/s)
        hand_stream_flowrate = 10/24.16/60 # SLPM to moles/s
        hand_stream = ct.importPhase('cantera_biomass/GasifierSpecies.cti')
        hand_stream.set(X = 'CO2:0.5, H2O:0.5', T = 400, P = 101325)
        tot_enth = hand_stream.enthalpy_mole()/1000/1000*hand_stream_flowrate
        
        self.assertAlmostEqual(test_stream.get_enthalpy('kJ/s'), tot_enth, 2)

class HeatExchangeTests(unittest.TestCase):
    """Needs to:
    0 Correctly calculate efficiency and effectiveness
    0 Raise exception if conversion, molar flowrates (inlet and outlet), or enthalpy calcs not done
    0 Raise exception on bad area or input power arguments
    """

    pass

class IntegralMeasurementTests(unittest.TestCase):
    pass

class UncertaintyTests(unittest.TestCase):
    """Needs to:
    0 Correctly add an uncertainty vector for a given calculation function given a set of uncertainties
    0 Raise an error if the uncertainty factors for the passed parameter list are not yet defined
    0 Raise an error on bad parameter types
    0 Raise an error if the function is not callable
    
    #FLESH THIS OUT LATER, ONCE EVERYTHING ELSE IS WORKING
    """
    """Carbon conversions must be calculated accurately and appropriate units set"""
    
    def testCorrectUncertaintyCalcs(self):
    
        gts = lfl.GasifierProcTS(start = datetime.datetime(1981,07,06,13,12,12),end=datetime.datetime(1981,07,06,13,12,16))
        for (key, value) in LoadDataTests.general_library.items():
            gts[key] = value
        gts.set_units(LoadDataTests.units)
        biomass_feed = lfl.Stream('ME_101',flowrate = gts.get_val('ME_101'),composition = {'biomass':1.00}, basis = "mass")
        entrainment_gas_feed = lfl.Stream('MFC_101', flowrate = gts.get_val('MFC_101'),composition = {'N2':0.95, 'Ar':0.05}, basis = "gas_volume")
        methane_gas_feed = lfl.Stream('MFC_102', flowrate = gts.get_val('MFC_102'), composition = {'CH4':1.00}, basis = "gas_volume")
        gas_exit = lfl.Stream('FE_101', flowrate = gts.get_val('FE_101'),basis = "gas_volume")
        
        biomass_feed.temperature = (25,"C")
        biomass_feed.pressure = (101325, "Pa")


        entrainment_gas_feed.temperature = (25,"C")
        methane_gas_feed.temperature = (25, "C")
        gas_exit.temperature = (25, "C")
        entrainment_gas_feed.pressure = (101325, "Pa")
        methane_gas_feed.pressure = (101325, "Pa")
        gas_exit.pressure = (101325, 'kg/s^2/m')

        exit_list = ['H2', 'CO', 'CO2', 'CH4', 'C2H6', 'N2', 'Ar']
        for specie in exit_list:
            key = "%s_GC" % specie
            gas_exit.composition[specie] = gts[key]
        gts.inlet_streams = [entrainment_gas_feed, methane_gas_feed, biomass_feed]
        gts.outlet_streams = [gas_exit]
        gts.proc_elements = ['C', 'H', 'O']
        gts.proc_species = ['H2', 'CO', 'CO2', 'CH4', 'C2H6', 'N2']
        biomass_breakdown = {}
        biomass_breakdown['biomass'] = {'C':0.5, 'O': 0.4, 'H': 0.1}
        biomass_feed.special_species = biomass_breakdown
        gts.generate_inlet_outlet_elemental_flows()
        gts.generate_inlet_outlet_species_flows()

        #Create uncertainty objects (like structs) to hold the uncertainty parameters for the critical numbers


#Unit tests added by Ryon 9/3/2013

class ProcessObjectTests(unittest.TestCase):
    """Needs to:
    0 Correctly calculate total inlet enthalpy.
    0 Correctly calculate total inlet entropy.
    0 Correctly calculate total outlet enthalpy.
    0 Correctly calculate total outlet entropy.
    0 Correctly calculate delta H.
    0 Correctly calculate delta S.
    
    0 Raise error if inlets are not in list object.
    0 Raise error if outlets are not in list object.
    0 Raise error if any inlet in inlets list not a Stream object.
    0 Raise error if any outlets in outlets list not a Stream object.
    """

    def setUp(self):
        pass
    def testInletEnthalpy(self):
        pass
    def testOutletEnthalpy(self):
        pass

class MixerTests(unittest.TestCase):
    """Needs to:
    0 Correctly calculate outlet composition of a mixer.
    0 Correctly calculate outlet temperature of a mixer.
    0 Correctly calculate outlet pressure of a mixer.
    """
    
    def setUp(self):
        pass
        

class SpaceTimeTests(unittest.TestCase):
    """Needs to:
    0 Correctly calculate space time
    """
    def setUp(self):
        gts = lfl.GasifierProcTS(start = datetime.datetime(1981,07,06,13,12,12), end = datetime.datetime(1981,07,06,13,12,16))
        for (key, value) in LoadDataTests.general_library.items():
            gts[key] = value
        gts.set_units(LoadDataTests.units)
        biomass_feed = lfl.Stream('biomass_feed', flowrate = gts.get_val('ME_101'), composition = {'biomass':1.00}, basis = "mass", temperature = (25, 'C'), pressure = (50, 'psig'))
        ent_1 = lfl.Stream('ent_1', flowrate = gts.get_val('MFC_102'), composition = {'N2':1.00}, basis = "std_gas_volume", temperature = (25, 'C'), pressure = (50, 'psig'))
        ent_2 = lfl.Stream('ent_2', flowrate = gts.get_val('MFC_103'), composition = {'CO2':1.00}, basis = "std_gas_volume", temperature = (25, 'C'), pressure = (50, 'psig'))
        ent_3 = lfl.Stream('ent_3', flowrate = gts.get_val('MFC_104'), composition = {'Ar':1.00}, basis = "std_gas_volume", temperature = (25, 'C'), pressure = (50, 'psig'))
        
        inlet_streams = [biomass_feed, ent_1, ent_2, ent_3]
        reactor_vol = (1.5*1.5*np.pi/4*24, 'in^3')
        
        # Hand calculate residence time...
        
        hand_vol = 1.5*1.5*np,pi/4*24*0.0163871 #Convert to liters
        

if __name__ == "__main__":
    
    unittest.main()
    
