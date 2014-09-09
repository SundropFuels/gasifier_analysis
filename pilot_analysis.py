#!/usr/bin/env python

from LabFunctionLib import *
import datetime
import numpy as np
import db_toolbox as db
import Thermo
import argparse
import csv
import getpass

class PilotDataAnalysis:
    """The basic data analysis class for pilot experiments"""

    def __init__(self, run_id, user, password, server = "localhost", run_information = None):
        #Create the pilot data frame, and load the data from the SQL database.


               
        self.interface_raw = db.db_interface(host = server, user = user, passwd = password)
        self.interface_raw.connect()
        q = db.use_Query("pilot_run_db")
        self.interface_raw.query(q)

        self.interface_proc = db.db_interface(host = server, user = user, passwd = password)
        self.interface_proc.connect()
        q = db.use_Query("pilot_proc_db")
        self.interface_proc.query(q)

        self.run_id = run_id

        self.run_info = RunInformation()
        #self.run_info.info = run_information
        self._load_run_info()
        self._load_timeseries_data()
        self._setup_standard_streams()
        
    def _load_run_info(self):

        #see if we can get a complete view from this
        self.run_info.SQL_load(self.interface_proc, table = 'run_plan_view', run_id = self.run_id)
        
        #set up the tube size
        if self.run_info.info['tube_diameter'] is not None:
            self.reactor_size = (self.run_info.info['tube_diameter']**2 * np.pi/4*12.0*18.0, 'in^3')
        else:
            self.reactor_size = (None, None)
        

    def _load_timeseries_data(self):
        """Loads the timeseries data into the database"""
        
        self.gts = GasifierProcTS(start = self.run_info.info['ss_start'], end = self.run_info.info['ss_stop'])
        self.gts.SQL_load(self.interface_raw,'analysis_view', glossary = 'glossary_tbl') #This line needs to automatically load the units
        #Need to build the glossary using the SQL tools
        q = db.select_Query(objects = ['tag_number', 'simple_name', 'units'], table = "glossary_tbl")
	glossary = self.interface_raw.query(q)
	self.glossary = {}
	self.gl_units = {}
	for row in glossary:
            self.glossary[row['tag_number']] = row['simple_name']
            self.gl_units[row['simple_name']] = row['units']
        
        self.gts.glossary_replace(self.glossary)
        
        self.gts.set_units(self.gl_units)
        self.gts.replace_None_w_nan_all()
        
    def _setup_standard_streams(self):
        """Sets up the standard streams and material compositions for analysis"""
        
        ########################################
        ### Following streams must be defined ##
        ########################################
        
        #1 Entrainment CO2
        #2 Sweep CO2
        #3 Steam
        #4 Argon
        #5 Biomass
        #6 Superheater N2 purge (technically, part of the steam, but adding separately for unit convenience)
        #7 Gas Exit
        
        #######
        ### ###
        #######
        
        #Setup the standard temperatures and pressures
        MFC_SP = (101325.0, 'Pa')
        MFC_ST = (70.0, 'F')
        
        #1 Entrainment        
        entrainment = Stream('entrainment', flowrate = self.gts.val_units('mass_flow_entrainment'), composition = {'CO2':1}, basis = 'std_gas_volume')
        entrainment.set_temperature((25.0, 'C')) #Assumed ambient
        entrainment.set_pressure(self.gts.val_units('pressure_bell_housing')) #PI_924502
        entrainment.std_temperature = MFC_ST
        entrainment.std_pressure = MFC_SP
        
        #2 Sweep
        sweep = Stream('sweep', flowrate = self.gts.val_units('mass_flow_sweep'), composition = {'CO2':1}, basis = 'std_gas_volume')
        sweep.set_temperature((25.0, 'C')) #Assumed ambient
        sweep.set_pressure(self.gts.val_units('pressure_bell_housing')) #PI_924502
        sweep.std_temperature = MFC_ST
        sweep.std_pressure = MFC_SP
        
        #3 Steam        
        steam = Stream('steam', flowrate = self.gts.val_units('mass_flow_steam'), composition = {'H2O':1}, basis = 'mass')
        steam.set_temperature(self.gts.val_units('temp_steam_gasifier_inlet')) #TI_983036
        steam.set_pressure(self.gts.val_units('pressure_bell_housing')) #
        
        #4 Argon
        argon_tracer_feed = Stream('argon_tracer', flowrate = self.gts.val_units('mass_flow_argon'), composition = {'Ar':1}, basis = 'std_gas_volume')
        argon_tracer_feed.set_temperature((25.0, 'C')) #Assumed ambient
        argon_tracer_feed.set_pressure(self.gts.val_units('pressure_bell_housing')) #PI_924502
        argon_tracer_feed.std_temperature = MFC_ST
        argon_tracer_feed.std_pressure = MFC_SP
                
        #5 Biomass
        ####Flowrate from the slope/intercept method

        biomass_flowrate = (self.gts['roto_feed_op'] * self.run_info['feeder_slope'] + self.run_info['feeder_intercept'], 'lb/hr')
        self.gts['biomass_flowrate'] = biomass_flowrate[0]
        self.gts.units['biomass_flowrate'] = biomass_flowrate[1]


        ####Flowrate from the mass loss of the feeder method
        fit = np.polyfit(self.gts['counter'], self.gts['mass_feed_vessel'], 1)  #Fit a line through the feed vessel data; slope is fit[0]
        biomass_flowrate_hopper_massloss = (fit[0]*-1.0*3600*np.ones(len(self.gts['mass_feed_vessel'])), 'lb/hr')
        self.gts['biomass_flowrate_hopper_massloss'] = biomass_flowrate_hopper_massloss[0]
        self.gts.units['biomass_flowrate_hopper_massloss'] = biomass_flowrate_hopper_massloss[1]

        #Alternately, we could get this flowrate from the feeder by looking only at where the flowrate drops
        lag1 = self.gts['mass_feed_vessel'].shift(-1)
        steps = lag1 - self.gts['mass_feed_vessel']
        self.gts['mass_feed_vessel_stepped'] = self.gts['mass_feed_vessel']
        self.gts['mass_feed_vessel_stepped'][steps<0.1] = np.nan  #Only take the places where it changes
        #now, we have two options -- we could try to do a derivative for instantaneous flowrate, or we could just fit a line
        #let's fit a line to start
        fit = np.polyfit(self.gts['counter'][np.isfinite(self.gts['mass_feed_vessel_stepped'])], self.gts['mass_feed_vessel_stepped'][np.isfinite(self.gts['mass_feed_vessel_stepped'])], 1)
        biomass_flowrate_hopper_massloss_s = (fit[0]*-1.0*3600*np.ones(len(self.gts['mass_feed_vessel_stepped'])), 'lb/hr')
        self.gts['biomass_flowrate_hopper_massloss_s'] = biomass_flowrate_hopper_massloss_s[0]
        self.gts.units['biomass_flowrate_hopper_massloss_s'] = biomass_flowrate_hopper_massloss_s[1]

 


        biomass_feed = Stream('biomass_feed',flowrate = biomass_flowrate_hopper_massloss_s, composition = {'H2O':self.run_info['moisture']/100.0, 'biomass':1.00-self.run_info['moisture']/100.0}, basis = "mass")
        biomass_feed.set_temperature((25.0, 'C'))
        biomass_feed.set_pressure(self.gts.val_units('pressure_bell_housing'))      

        #6 N2 Superheater Purge
        N2_superheater_purge = Stream('N2_superheater_purge', flowrate = self.gts.val_units('mass_flow_superheater_purge'), composition = {'N2':1.0}, basis = "std_gas_volume")
        N2_superheater_purge.set_temperature(self.gts.val_units('temp_steam_gasifier_inlet'))
        N2_superheater_purge.set_pressure(self.gts.val_units('pressure_bell_housing'))
        N2_superheater_purge.set_std_temperature = MFC_ST
        N2_superheater_purge.set_std_pressure = MFC_SP

        #Set up exit gas flowrates
        gas_exit = self.gts.outlet_stream_from_tracer([argon_tracer_feed],"Molar", "Ar", self.gts['Ar_MS']/100.0, 'gas_exit')
        
        #Need to set up the exit gas composition now -- we can actually build a variety of different streams here and use them as necessary 
        
        mass_spec_list = ['N2','Ar','H2O','H2','CO','CO2','CH4','C2H6','C2H4','C2H2','C3H8','C3H6','C4H8','C4H10','CH3CHCH3CH3','C6H6','C7H8','C6H4CH3CH3','C6H5CH2CH3','C10H8','H2S']
        
        composition = {}
        for specie in mass_spec_list:
            composition[specie] = self.gts['%s_MS' % specie]/100.0
        #!!#
        ppm_list = ['C6H6','H2S','C10H8','C7H8','C6H4CH3CH3','C6H5CH2CH3']
        for key in ppm_list:
            composition[key] /= 10000.0
        #!!#
        gas_exit.set_pressure(self.gts.val_units('pressure_outlet'))
        gas_exit.set_composition(composition)
        gas_exit.set_temperature(self.gts.val_units('temp_gasifier_exit'))
        self.gts['exit_gas_flowrate'] = gas_exit.flowrate[0]
        self.gts.units['exit_gas_flowrate'] = 'mol/s'
               

        #!!# This is where you will include the laser Raman data, if desired


        self.gts.inlet_streams = [entrainment, sweep, steam, argon_tracer_feed, biomass_feed, N2_superheater_purge]
        
        self.gts.outlet_streams = [gas_exit]
        
        self.gts.proc_elements = ['C', 'H', 'O']
        self.gts.proc_species = ['N2','Ar','H2O','H2','CO','CO2','CH4','C2H6','C2H4','C2H2','C3H8','C3H6','C4H8','C4H10','CH3CHCH3CH3','C6H6','C7H8','C6H4CH3CH3','C6H5CH2CH3','C10H8','H2S']
        self.gts.inert_species = ['N2', 'Ar']
        
        #set up the biomass information from the run table
        biomass_breakdown = {}
        biomass_breakdown['biomass'] = {}
        for item in ['c', 'h', 'n']:
            biomass_breakdown['biomass'][item.upper()] = self.run_info['w_%s' % item]/100.0
        
        #Following line overrides assignment for individual sample carbon content and uses bag average carbon content instead. 
        #biomass_breakdown['biomass']['C'] = self.run_info['bag_c']/100
        
        #End additional row
        biomass_breakdown['biomass']['O'] = 1 - sum(biomass_breakdown['biomass'].values())
        biomass_feed.special_species = biomass_breakdown
        
        

    def calculate_standard_things(self):
        """Runs the standard calculations on the data"""
        #1. Calculate the elemental and species flows
        self.gts.generate_inlet_outlet_elemental_flows()
        self.gts.generate_inlet_outlet_species_flows()
        print "species flows calculated"
        #2. Calculate carbon conversions
        self.gts.generate_carbon_conversions()
        self.gts.generate_C_mass_balance()
        self.gts.generate_CH4_yield()
        print "carbon conversion calculated"
        #3. Calculate changes in enthalpy and entropy
        self.gts.calc_max_dH(temperature = [self.gts['temp_skin_tube_middle'],self.gts.units['temp_skin_tube_middle']], pressure = [self.gts['pressure_outlet'],self.gts.units['pressure_outlet']], units = 'kW')
        self.gts.generate_enthalpy_change('kW')
        #self.gts.generate_entropy_change(self, 'kW/K')

        #4. Normalize the compositions
        self.gts.generate_normalized_compositions()
        print "compositions normalized"
        #calculate inlet partial pressures
        self.gts['pp_CO2'] = self.gts['CO2_inlet']/(self.gts['CO2_inlet']+self.gts['Ar_inlet']+self.gts['H2O_inlet']+self.gts['N2_inlet'])*(self.gts['pressure_ako']+14.7)
        self.gts['pp_H2O'] = self.gts['H2O_inlet']/(self.gts['CO2_inlet']+self.gts['Ar_inlet']+self.gts['H2O_inlet']+self.gts['N2_inlet'])*(self.gts['pressure_ako']+14.7)
        self.gts['pp_Ar'] = self.gts['Ar_inlet']/(self.gts['CO2_inlet']+self.gts['Ar_inlet']+self.gts['H2O_inlet']+self.gts['N2_inlet'])*(self.gts['pressure_ako']+14.7)
        
        #5. Calculate tar loads
        self.gts.calc_tar_rate(self.gts.outlet_streams[0], tar_list = ['C6H6', 'C7H8', 'C10H8', 'C6H4CH3CH3', 'C6H5CH2CH3'], inclusive_tar_list = ['C2H2', 'C2H4', 'C2H6', 'C3H8', 'C3H6', 'C4H8', 'C4H10', 'C6H6', 'C7H8', 'C10H8', 'CH3CHCH3CH3', 'C6H4CH3CH3', 'C6H5CH2CH3'])
        print "tar loads calculated"
        #6. Calculate the space time
        self.gts.calc_space_time(self.reactor_size, 'biomass')
        print "space time calculated"
        #7. Calculate optical thickness

        sizes = ['10','50','90']
        for size in sizes:
            if self.run_info.info['d%s'%size] is None:
                self.run_info.info['d%s'%size] = np.nan

        tube_id = (self.run_info['tube_diameter'], 'in')

        self.gts.calc_optical_thickness(tubeD = tube_id, density = (1400, 'kg/m^3'), 
                                        particle_size = {'d10':(self.run_info.info['d10']*10**-6, 'm'), 
                                                         'd50':(self.run_info.info['d50']*10**-6, 'm'), 
                                                         'd90':(self.run_info.info['d90']*10**-6, 'm')})

        print "optical thickness calculated"
        #8. Calculate integral measures
        self.gts.generate_averages_and_stdevs(cols = self.get_analysis_config_list())
        
        #!!# New stuff

        #9. Calculate uncertainties
        ####NOT IMPLEMENTED YET
       
    def generate_output_file(self, filename):
        """Writes a csv file to filename from the timeseries data"""
        #Meta-data first
        f = open(filename, 'w')
        f.write('Start time: %s\n' % self.run_info['ss_start'])  #May need to do some formatting here
        f.write('End time: %s\n' % self.run_info['ss_stop']) #May need to do some formatting here
        #f.write('Biomass type', self.run_info.biomass_type)
        #f.write('Biomass moisture', self.run_info.moisture)

        #NEED TO INCLUDE SETPOINT INFORMATION HERE, POSSIBLY

        #now, write the data to the csv file
        f.close()
        #self.gts.write_csv(filename, mode = 'append', col_list = ["timestamp", "X_tot", "X_good", "X_std", "CO2_normalized", "CO_normalized", "CH4_normalized", "H2_normalized", "C2H2_outlet", "C2H4_outlet", "CH4_outlet", "C6H6_outlet", "tar_loading", "tar_loading_incl"])
        self.gts.write_csv(filename, mode = 'append')

    def upload_to_database(self):

        #upload the time series data
        
        self.gts.SQL_db_upload(self.interface_proc, table = "pilot_proc_data_tbl")

        print "completed upload of timeseries data"
        #upload the integral data
        self.upload_integral_data()

    def upload_integral_data(self):
        #need to build the query elements
        objects = {}
        objects['run_id'] = str(self.run_id)
        objects['N_total'] = str(len(self.gts.index))
        objects['N_MS'] = str(len(self.gts['CO_MS'][np.isfinite(self.gts['CO_MS'].astype(float))]))
        objects['analysis_ts'] = datetime.datetime.strftime(datetime.datetime.today(), '%Y-%m-%d %H:%M:%S')
        for key in self.gts.avgs:
            objects['%s_avg' % key] = str(self.gts.avgs[key])
        for key in self.gts.stdevs:
            objects['%s_std' % key] = str(self.gts.stdevs[key])
        
        #try to upload the data
        try:
            q = db.insert_Query(objects, table = "integral_tbl")
            #print q.getQuery()
            self.interface_proc.query(q)
        except db.DBToolboxError:
            try:
                
                q = db.update_Query(objects, table = "integral_tbl", conditions = ["run_id = %s" % objects['run_id']])
                #print q.getQuery()
                self.interface_proc.query(q)
            except db.DBToolboxError:
                print "Crashed trying to upload to integral data table"

    def get_analysis_config_list(self):
        #q = db.use_Query("lab_proc_db")
        objects = ['avg_std_cols']
        q = db.select_Query(objects, table = 'analysis_config_tbl', condition_list = ["active = 1"])
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
    parser.add_argument('--user', type=str,action='store')
    parser.add_argument('--pswd', type=str, action='store')
    parser.add_argument('--host', type=str, action ='store')
    args = parser.parse_args()

    if args.user is not None:
        user = args.user
    if args.pswd is not None:
        pswd = args.pswd
    if args.host is not None:
        host = args.host
    else:
        host = 'localhost'


    if args.user is None or args.pswd is None:
        user = raw_input('User: ')
        pswd = getpass.getpass()
        
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
        analyzer = PilotDataAnalysis(run_id = run_id, user = user, password = pswd, server = host)
        print "Data loaded"
        analyzer.calculate_standard_things()
        print "Standard things calculated"
        #analyzer.generate_output_file('run100.csv')
        analyzer.upload_to_database()
        print "Data uploaded to database"
