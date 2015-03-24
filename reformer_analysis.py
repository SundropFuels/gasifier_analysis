#!/usr/bin/env python

from LabFunctionLib import *
import datetime
import numpy as np
import db_toolbox as db
import Thermo
import argparse
import csv
import getpass
import unitConversion as uc
conv = uc.UnitConverter()

### This is a newer file than gasifier_analysis.py.  Because of time constraints on these experiments, data upload to databases was not set up fully.  However, there are some more elegant solutions in this file compared to gasifier_analysis.py.

class GasifierDataAnalysis:
    """The basic data analysis class for gasifier experiments"""

    def __init__(self, run_id, user, password, run_information = None):
        #Create the gasifier data frame, and load the data from the SQL database (this will be automated through the interface later)
               
        self.interface_raw = db.db_interface(host = "192.168.13.51", user = user, passwd = password)
        self.interface_raw.connect()
        q = db.use_Query("lab_run_db")
        self.interface_raw.query(q)

        self.interface_proc = db.db_interface(host = "192.168.13.51", user = user, passwd = password)
        self.interface_proc.connect()
        q = db.use_Query("lab_proc_db")
        self.interface_proc.query(q)

        self.run_id = run_id

        self.run_info = RunInformation()
        #self.run_info.info = run_information
        self._load_run_info()
        self._load_timeseries_data()
        self._setup_standard_streams()
        
    def _load_run_info(self):

        #see if we can get a complete view from this
        self.run_info.SQL_load(self.interface_proc, table = 'rfmr_run_plan_view', run_id = self.run_id)
        #set up the tube size
        self.reactor_size = (24.0*self.run_info.info['tube_in_dia']**2*np.pi/4, 'in^3')
        self.tube_id = (self.run_info.info['tube_in_dia'], 'in')
        self.reactor_sa = self.tube_id[0]*np.pi*24*0.00064516
        

    def _load_timeseries_data(self):
        """Loads the timeseries data into the database"""
        
        self.gts = GasifierProcTS(start = self.run_info.info['ss_start'], end = self.run_info.info['ss_stop'])
        self.gts.SQL_load(self.interface_raw,'rfmr_analysis_view') #This line needs to automatically load the units
        #Need to build the glossary using the SQL tools
        q = db.select_Query(objects = ['tag_name', 'simple_name', 'units'], table = "rfmr_glossary_tbl")
	glossary = self.interface_raw.query(q)
	self.glossary = {}
	self.gl_units = {}
	for row in glossary:
            self.glossary[row['tag_name']] = row['simple_name']
            self.gl_units[row['simple_name']] = row['units']

        self.gts.glossary_replace(self.glossary)

        self.gts.set_units(self.gl_units)
        self.gts.replace_None_w_nan_all()

        # Filter out mass flow controller PV's which give false or negative feed rates with zero flow.  This is a better solution compared to gasifier_analysis.py, as it accounts for possibility of an operator mis-setting the setpoint on the HMI.
        for i in [i for i in self.gts.columns if i.startswith('mass_flow_')]:
            self.gts[i][self.gts[i] < 0.5] = 0
        # Get list of columns in raw data frame for uploading calculated columns later
        self.raw_cols = self.gts.columns

    def _setup_standard_streams(self):
        """Sets up the standard streams and material compositions for analysis"""

        conv = uc.UnitConverter()
        # Establish MFC_101 gas type
        if self.gts['MFC_101_GAS'][0] == 0: 
            gas_101 = "CO2"
        elif self.gts['MFC_101_GAS'][0] == 1:
            gas_101 = "N2"
        else:
            gas_101 = "N2"
        # Establish MFC_201 gas type    
        if self.gts['MFC_201_GAS'][0] == 0: 
            gas_201 = "CO2"
        elif self.gts['MFC_201_GAS'][0] == 1:
            gas_201 = "N2"
        else: 
            gas_201 = "CO2"
        # Apply MFC_101 correction factor y=ax+b, default is a=1, b=0.  Accounts for mis-calibration issues.
        a_101 = self.run_info.info['MFC_101_slope']
        b_101 = self.run_info.info['MFC_101_int']
        self.gts['mass_flow_101'] = a_101 * self.gts['mass_flow_101'] + b_101
        # Apply MFC_201 correction factor y=ax+b, default is a=1, b=0 .  Accounts for mis-calibration issues.
        a_201 = self.run_info.info['MFC_201_slope']
        b_201 = self.run_info.info['MFC_201_int']
        self.gts['mass_flow_201'] = a_201 * self.gts['mass_flow_201'] + b_201

	# Set up gas feed streams.  These streams are set up cleaner than gasifier_analysis.py 
        # Inert gas feed streams
        feed_101 = Stream('101_feed', flowrate = self.gts.val_units('mass_flow_101'), composition = {gas_101:1.00}, temperature = (25, 'C'), pressure = self.gts.val_units('pressure_gas_inlet'), basis = 'std_gas_volume')
        Ar_feed = Stream('Ar_feed', flowrate = self.gts.val_units('mass_flow_Ar'), composition = {'Ar':1.00}, temperature = (25, 'C'), pressure = self.gts.val_units('pressure_gas_inlet'), basis = "std_gas_volume")
        feed_201 = Stream('201_feed', flowrate = self.gts.val_units('mass_flow_201'),composition = {gas_201:1.00}, temperature = (25, 'C'), pressure = self.gts.val_units('pressure_gas_inlet'), basis = "std_gas_volume")

        # Combustible gas feed streams
        CO_feed = Stream('CO_feed', flowrate = self.gts.val_units('mass_flow_CO'), composition = {'CO':1.00}, temperature = (25, 'C'), pressure = self.gts.val_units('pressure_flammable_inlet'), basis = 'std_gas_volume')
        CH4_feed = Stream('CH4_feed', flowrate = self.gts.val_units('mass_flow_CH4'), composition = {'CH4':1.00}, temperature = (25, 'C'), pressure = self.gts.val_units('pressure_gas_inlet'), basis = 'std_gas_volume')
        C2H2_feed = Stream('C2H2_feed', flowrate = self.gts.val_units('mass_flow_C2H2'), composition = {'C2H2':1.00}, temperature = (25, 'C'), pressure = self.gts.val_units('pressure_gas_inlet'), basis = 'std_gas_volume')
        C2H4_feed = Stream('C2H4_feed', flowrate = self.gts.val_units('mass_flow_C2H4'), composition = {'C2H4':1.00}, temperature = (25, 'C'), pressure = self.gts.val_units('pressure_gas_inlet'), basis = 'std_gas_volume')
        C2H6_feed = Stream('C2H6_feed', flowrate = self.gts.val_units('mass_flow_C2H6'), composition = {'C2H6':1.00}, temperature = (25, 'C'), pressure = self.gts.val_units('pressure_gas_inlet'), basis = 'std_gas_volume')

        # Combustible liquid feed streams, convert to mass flow from volumetric flow
#        self.gts.convert_col_units('volumetric_flow_benzene_sp', 'mL/min')
#        self.gts['C6H6_feed'] = self.gts['volumetric_flow_benzene_sp'] * 0.87650
#        self.gts.units['C6H6_feed'] = 'g/min'
#        C6H6_feed = Stream('C6H6_feed', flowrate = self.gts.val_units('C6H6_feed'), composition = {'C6H6':1}, temperature = self.gts.val_units('temperature_vapor_to_lance_inlet'), pressure = self.gts.val_units(''), basis = 'mass')

#        self.gts.convert_col_units('volumetric_flow_benzene_sp', 'mL/min')
#        self.gts['C10H8_feed'] = self.gts['volumetric_flow_naphthalene_sp'] * 1.14
#        self.gts.units['C10H8_feed'] = 'g/min'
#        C6H6_feed = Stream('C10H8_feed', flowrate = self.gts.val_units('C10H8_feed'), composition = {'C10H8':1}, temperature = self.gts.val_units('temperature_vapor_to_lance_inlet'), pressure = self.gts.val_units(''), basis = 'mass')

        # Steam stream, convert steam to mass flow from volumetric flow
        self.gts.convert_col_units('volumetric_flow_steam_water_sp', 'mL/hr')
        self.gts['H2O_feed'] = self.gts['volumetric_flow_steam_water_sp']
        self.gts.units['H2O_feed'] = 'g/hr'
        self.gts.convert_col_units('H2O_feed', 'lb/hr')
        steam_feed = Stream('steam_feed', flowrate = self.gts.val_units('H2O_feed'), composition = {'H2O':1.00}, temperature = self.gts.val_units('temperature_vapor_to_lance_inlet'), pressure = self.gts.val_units('pressure_steam_water_supply'), basis = 'mass')

        MFC_SP = (101325.0, 'Pa')
        MFC_ST = (70.0, 'F')

        #Setup the standard temperatures and pressures
        feed_101.std_temperature = MFC_ST
        Ar_feed.std_temperature = MFC_ST
        CO_feed.std_temperature = MFC_ST
        feed_201.std_temperature = MFC_ST
        CH4_feed.std_temperature = MFC_ST
        C2H2_feed.std_temperature = MFC_ST
        C2H4_feed.std_temperature = MFC_ST
        C2H6_feed.std_temperature = MFC_ST

        feed_101.std_pressure = MFC_SP
        Ar_feed.std_pressure = MFC_SP
        CO_feed.std_pressure = MFC_SP
        feed_201.std_pressure = MFC_SP
        CH4_feed.std_pressure = MFC_SP
        C2H2_feed.std_pressure = MFC_SP
        C2H4_feed.std_pressure = MFC_SP
        C2H6_feed.std_pressure = MFC_SP

        #Set up exit gas flowrates

        gas_exit = self.gts.outlet_stream_from_tracer([Ar_feed],"Molar", "Ar", self.gts['Ar_MS']/100.0, 'gas_exit')
        #Need to set up the exit gas composition now -- we can actually build a variety of different streams here and use them as necessary
        self.mass_spec_list = ['N2', 'Ar', 'H2', 'CO', 'CO2', 'CH4', 'C2H6', 'C2H4', 'C2H2', 'C3H8', 'C3H6', 'C6H6', 'C7H8', 'C10H8', 'H2O']

        composition = {}
        for specie in self.mass_spec_list:
            composition[specie] = conv.convert_units(self.gts['%s_MS' % specie], self.gts.units['%s_MS' % specie], 'fraction')
#        Will not need following lines if above unit conversion works as planned
#        ppm_list = ['C6H6','H2S','C10H8','C7H8']
#        for key in ppm_list:
#            composition[key] /= 10000.0
        gas_exit.set_pressure(self.gts.val_units('pressure_ash_knockout'))
        gas_exit.set_composition(composition)
        gas_exit.set_temperature(self.gts.val_units('temperature_reactor_skin_middle'))
        self.gts['exit_gas_flowrate'] = gas_exit.flowrate[0]
        self.gts.units['exit_gas_flowrate'] = 'mol/s'

        #interpolate the Ar column to get measurements for GC comparisons
#        GC_list = ['H2', 'CO', 'CO2', 'CH4', 'Ar', 'N2']
#        for specie in GC_list:
#            self.gts['%s_GC_interp' % specie] = self.gts.interpolate_col('counter', '%s_GC' % specie)

#        gas_exit_GC = self.gts.outlet_stream_from_tracer([Ar_feed], 'Molar', 'Ar', self.gts['Ar_GC_interp']/100.0, 'gas_exit_interp')
#        gas_exit_GC.set_pressure(self.gts.val_units('pressure_product_gas_downstream_filters'))
#        gas_exit_GC.set_temperature(self.gts.val_units('temp_exit_gas'))
#        gc_composition = {}
#        for specie in GC_list:
#            gc_composition[specie] = self.gts["%s_GC_interp" % specie]
#            gas_exit_GC.set_composition(gc_composition)

        self.gts.inlet_streams = [Ar_feed, feed_101, feed_201, CH4_feed, steam_feed]
        
        self.gts.outlet_streams = [gas_exit]
        
        self.gts.proc_elements = ['C', 'H', 'O']
        self.gts.proc_species = ['H2', 'CO', 'CO2', 'CH4', 'C2H6', 'N2', 'C2H4', 'C2H2', 'C6H6', 'Ar', 'C7H8', 'C10H8', 'C3H8', 'C3H6','H2O']
        self.gts.inert_species = ['N2', 'Ar']

    def calculate_standard_things(self):
        """Runs the standard calculations on the data"""
        #1. Calculate the elemental and species flows
        self.gts.generate_inlet_outlet_elemental_flows()
        self.gts.generate_inlet_outlet_species_flows()

        #2. Calculate carbon conversions and mass balances.  These had to be defined manually in this file because there was little time to create them properly in labfunctionlib.  
        self.gts['X_CH4'] = 1-self.gts['CH4_outlet']/self.gts['CH4_inlet']
        self.gts['X_C6H6'] = 1 - self.gts['C6H6_outlet']/self.gts['C6H6_outlet']
        self.gts['Y_CO_CO2'] = (self.gts['CO_outlet'] + self.gts['CO2_outlet'] - self.gts['CO2_inlet'])/self.gts['CH4_inlet']
        self.gts['Y_soot'] = (self.gts['C_inlet'] - self.gts['C_outlet'])/self.gts['CH4_inlet']
        self.gts['Y_C6H6'] = self.gts['C6H6_outlet']/self.gts['CH4_inlet']
        self.gts['Y_C7H8'] = self.gts['C7H8_outlet']/self.gts['CH4_inlet']
        self.gts['Y_C10H8'] = self.gts['C10H8_outlet']/self.gts['CH4_inlet']
        self.gts['H2O:CH4'] = self.gts['H2O_inlet']/self.gts['CH4_inlet']

        #3. Calculate changes in enthalpy
        self.gts.calc_max_dH(temperature = [self.gts['temperature_reactor_skin_middle'],self.gts.units['temperature_reactor_skin_middle']], pressure = [self.gts['pressure_ash_knockout'],self.gts.units['pressure_ash_knockout']], units = 'kW')
        self.gts['dH/A'] = self.gts['dH_max']/self.reactor_sa
#        self.gts.generate_enthalpy_change('kW')

        #4. Normalize the compositions
        self.gts.generate_normalized_compositions()

        #calculate partial pressures

        total_inlet_flow = 0
        total_outlet_flow = 0
        for specie in self.mass_spec_list:
            total_inlet_flow += self.gts['%s_inlet' %specie]
            total_outlet_flow += self.gts['%s_outlet' %specie]

        for specie in self.mass_spec_list:
            self.gts['pp_%s_inlet' %specie] = self.gts['%s_inlet' %specie] / total_inlet_flow * (self.gts['pressure_ash_knockout'] + 14.7)
            self.gts['pp_%s_outlet' %specie] = self.gts['%s_outlet' %specie] / total_outlet_flow * (self.gts['pressure_ash_knockout'] + 14.7)

        #5. Calculate tar loads
        self.gts.calc_tar_rate(self.gts.outlet_streams[0])

        #6. Calculate the space time
        self.gts.calc_space_time(self.reactor_size, excluded_species = [None])

        #7. Calculate integral measures
#        self.gts.generate_averages_and_stdevs(cols = self.get_analysis_config_list())

    def generate_output_file(self, filename):
        """Writes a csv file to filename from the timeseries data"""
        #Meta-data first
        f = open(filename, 'w')
        f.write('Start time: %s\n' % self.run_info['ss_start']) #May need to do some formatting here
        f.write('End time: %s\n' % self.run_info['ss_stop']) #May need to do some formatting here

        #NEED TO INCLUDE SETPOINT INFORMATION HERE, POSSIBLY

        #now, write the data to the csv file
        f.close()
        #self.gts.write_csv(filename, mode = 'append', col_list = ["timestamp", "X_tot", "X_good", "X_std", "CO2_normalized", "CO_normalized", "CH4_normalized", "H2_normalized", "C2H2_outlet", "C2H4_outlet", "CH4_outlet", "C6H6_outlet", "tar_loading", "tar_loading_incl"])
        self.gts.write_csv(filename, mode = 'append')

    def upload_to_database(self):

        #upload the time series data
        self.gts.SQL_db_upload(self.interface_proc, table = "rfmr_proc_data_tbl")

        print "completed upload of timeseries data"
        #upload the integral data
        self.upload_integral_data()

    def upload_integral_data(self):
        #need to build the query elements
        objects = {}
        objects['run_id'] = str(self.run_id)
        objects['N_total'] = str(len(self.gts.index))
        objects['N_MS'] = str(len(self.gts['CO_MS'][np.isfinite(self.gts['CO_MS'].astype(float))]))
        objects['analysis_timestamp'] = datetime.datetime.strftime(datetime.datetime.today(), '%Y-%m-%d %H:%M:%S')
        for key in self.gts.avgs:
            objects['%s_avg' % key] = str(self.gts.avgs[key])
        for key in self.gts.stdevs:
            objects['%s_std' % key] = str(self.gts.stdevs[key])

        #try to upload the data
        try:
            q = db.insert_Query(objects, table = "rfmr_integral_tbl")
            #print q.getQuery()
            self.interface_proc.query(q)
        except db.DBToolboxError:
            try:

                q = db.update_Query(objects, table = "rfmr_integral_tbl", conditions = ["run_id = %s" % objects['run_id']])
                #print q.getQuery()
                self.interface_proc.query(q)
            except db.DBToolboxError:
                print "Crashed trying to upload to integral data table"

    def get_analysis_config_list(self):
        #q = db.use_Query("lab_proc_db")
        objects = ['avg_std_cols']
        q = db.select_Query(objects, table = 'rfmr_analysis_config_tbl', condition_list = ["active = 1"])
        col_list = []
        col_tuple = self.interface_proc.query(q)
        for d in col_tuple:
            col_list.append(d['avg_std_cols'])

        return col_list

    def generate_standard_report(self, rpt_fmt):
        #generate a standard report -- maybe given an xml format file to direct what goes into the report
        pass

def parse_list(txt):
    main_list = txt.split(",")
    run_id_list = []
    for sublist in main_list:
        if ":" in sublist:
            left = sublist.split(":")[0]
            right = sublist.split(":")[1]
            run_id_list.extend(range(int(left), int(right)+1))
        else:
            run_id_list.append(int(sublist))

    return run_id_list  

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = "Run a gasifier analysis")
    parser.add_argument('--run_range', type=str, action = 'store')
    parser.add_argument('--run_id',type=int, action ='store')
    parser.add_argument('--file',type=str,action = 'store')
    parser.add_argument('--csv', type = str, action = 'store')
    args = parser.parse_args()

    user = raw_input('User: ')
    pswd = getpass.getpass()

    results_df = pd.DataFrame()

    if args.run_id is not None:
        run_id_list = [args.run_id]

    elif args.run_range is not None:
        run_id_list = parse_list(args.run_range)

    elif args.file is not None:
        f = open(args.file)
        a = f.readline()
        a = a[:-1]

        run_id_list = parse_list(a)

    for run_id in run_id_list:
        print "Analyzing run %s..." % run_id

        analyzer = GasifierDataAnalysis(run_id = run_id, user = user, password = pswd)
        analyzer.calculate_standard_things()
        run_df = analyzer.gts.mean()
        results_df = results_df.append(run_df, ignore_index = True)
        print "Data loaded"
        print "Standard things calculated"
    results_df.index = run_id_list
    results_df.to_csv(args.csv)
#        analyzer.generate_output_file('run100.csv')
#        analyzer.upload_to_database()
#        print "Data uploaded to database for run_id %s" %run_id

