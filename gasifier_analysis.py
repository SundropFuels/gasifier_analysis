#!/usr/bin/env python

from LabFunctionLib import *
import datetime
import numpy as np
import db_toolbox as db
import Thermo
import argparse
import csv

class GasifierDataAnalysis:
    """The basic data analysis class for gasifier experiments"""

    def __init__(self, run_id, run_information = None):
        #Create the gasifier data frame, and load the data from the SQL database (this will be automated through the interface later)
        self.interface_raw = db.db_interface(host = "192.168.13.15", user = "chris", passwd = "udflyer87")
        self.interface_raw.connect()
        q = db.use_Query("lab_run_db")
        self.interface_raw.query(q)

        self.interface_proc = db.db_interface(host = "192.168.10.20", user = "chris", passwd = "cmp87ud01")
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
        self.run_info.SQL_load(self.interface_proc, table = 'gas_run_info_tbl', run_id = self.run_id)
        #set up the tube size
        if self.run_info.info['tube_dia'] == 2.0:
            self.reactor_size = (24.0*1.5, 'in^3')
        elif self.run_info.info['tube_dia'] == 1.5:
            self.reactor_size = (24*1.0, 'in^3')
        else:
            self.reactor_size = (None, None)
        
        

        #WILL NEED TO DO THE APPROPRIATE GLOSSARY SWITCHES OR JUST RENAME THINGS IN RUN_INFO TO THE RIGHT STUFF

    def _load_timeseries_data(self):
        """Loads the timeseries data into the database"""
        
        self.gts = GasifierProcTS(start = self.run_info.info['ss_start'], end = self.run_info.info['ss_stop'])
        self.gts.SQL_load(self.interface_raw,'gasifier_lv_GC_view') #This line needs to automatically load the units
        #Need to build the glossary using the SQL tools
        q = db.select_Query(objects = ['tag_number', 'simple_name', 'units'], table = "tag_glossary_tbl")
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
        biomass_feed = Stream('biomass_feed',flowrate = self.gts.val_units('mass_flow_brush_feeder'),composition = {'H2O':self.run_info['moisture']/100.0, 'biomass':1.00-self.run_info['moisture']/100.0}, basis = "mass")
        
        
        if self.gts['entrainment_gas_type'][0] == 0:
            e_type = "N2"
        elif self.gts['entrainment_gas_type'][0] == 1:
            e_type = "CO2"
        

        if self.run_info.info['downbed_gas_type'] == 0:
            db_type = "N2"
        elif self.run_info.info['downbed_gas_type'] == 1:
            db_type = "CO2"
        else:
            db_type = "CO2"


        """
        if self.gts['makeup_gas_type'][0] == 0:
            m_type = "N2"
        elif self.gts['makeup_gas_type'][0] == 1:
            m_type = "CO2"
        """
        m_type = "CO2"
        ###this is a problem!!!!###
        cross_brush_feed = Stream('cross_brush_feed', flowrate = self.gts.val_units('mass_flow_entrainment'), composition = {e_type:1.0}, basis = "std_gas_volume")
        if self.run_info.info['downbed_flow_rate'] < 0.01:
            down_bed_feed = Stream('down_bed_feed', flowrate = (self.gts['mass_flow_down_brush']*0.0,'L/min'), composition = {db_type:1.0}, basis = "std_gas_volume")
        else:
            down_bed_feed = Stream('down_bed_feed', flowrate = self.gts.val_units('mass_flow_down_brush'), composition = {db_type:1.0}, basis = "std_gas_volume")

        if self.run_info.info['makeup_flow'] < 0.01:
            makeup_feed = Stream('makeup_feed', flowrate = (self.gts['mass_flow_feed_vessel_pressure']*0.0,'L/min'), composition = {m_type:1.0}, basis = "std_gas_volume")
        else:
            makeup_feed = Stream('makeup_feed', flowrate = self.gts.val_units('mass_flow_feed_vessel_pressure'), composition = {m_type:1.0}, basis = "std_gas_volume")

        argon_tracer_feed = Stream('argon_tracer', flowrate = self.gts.val_units('mass_flow_argon_tracer'), composition = {'Ar':1.0}, basis = "std_gas_volume")
        methane_gas_feed = Stream('methane_feed', flowrate = self.gts.val_units('mass_flow_methane'),composition = {'CH4':1.00}, basis = "std_gas_volume")
        

        


        
        #Convert steam to mass flow 
        self.gts.convert_col_units('setpoint_steam_HPLC_pump', 'mL/hr')
        self.gts['steam_flow'] = self.gts['setpoint_steam_HPLC_pump']
        self.gts.units['steam_flow'] = 'g/hr'
        self.gts.convert_col_units('steam_flow', 'lb/hr')
        steam_feed = Stream('steam_feed', flowrate = self.gts.val_units('steam_flow'), composition = {'H2O':1.00}, basis = "mass")
        steam_feed.temperature = (self.gts['temp_steam_reactor_entry'], self.gts.units['temp_steam_reactor_entry'])
        steam_feed.pressure = (self.gts['pressure_reactor_gas_inlet'], self.gts.units['pressure_reactor_gas_inlet'])        
        
        MFC_SP = (101325.0, 'Pa')
        MFC_ST = (70.0, 'F')

        biomass_feed.temperature = (25.0, 'C')
        cross_brush_feed.temperature = (25.0, 'C')
        down_bed_feed.temperature = (25.0, 'C')
        makeup_feed.temperature = (25.0, 'C')
        argon_tracer_feed.temperature = (25.0, 'C')
        methane_gas_feed.temperature = (25.0, 'C')
        
        
        steam_feed.temperature = self.gts.val_units('temp_steam_reactor_entry')
        
        biomass_feed.pressure = self.gts.val_units('pressure_reactor_gas_inlet')
        cross_brush_feed.pressure = self.gts.val_units('pressure_reactor_gas_inlet')
        down_bed_feed.pressure = self.gts.val_units('pressure_reactor_gas_inlet')
        makeup_feed.pressure = self.gts.val_units('pressure_reactor_gas_inlet')
        argon_tracer_feed.pressure = self.gts.val_units('pressure_reactor_gas_inlet')
        methane_gas_feed.pressure = self.gts.val_units('pressure_reactor_gas_inlet')
        
        steam_feed.pressure = self.gts.val_units('pressure_reactor_gas_inlet')

        #Setup the standard temperatures and pressures
        cross_brush_feed.std_temperature = MFC_ST
        down_bed_feed.std_temperature = MFC_ST
        makeup_feed.std_temperature = MFC_ST
        argon_tracer_feed.std_temperature = MFC_ST
        methane_gas_feed.std_temperature = MFC_ST


        cross_brush_feed.std_pressure = MFC_SP
        down_bed_feed.std_pressure = MFC_SP
        makeup_feed.std_pressure = MFC_SP
        argon_tracer_feed.std_pressure = MFC_SP
        methane_gas_feed.std_pressure = MFC_SP

        #Set up exit gas flowrates
        
        
        gas_exit = self.gts.outlet_stream_from_tracer([argon_tracer_feed],"Molar", "Ar", self.gts['Ar_MS']/100.0, 'gas_exit')
        #Need to set up the exit gas composition now -- we can actually build a variety of different streams here and use them as necessary 
        mass_spec_list = ['C2H2', 'Ar', 'C6H6', 'CO2', 'C2H6', 'C2H4', 'H2S', 'H2', 'CH4', 'C10H8', 'N2', 'C3H8', 'C3H6', 'C7H8', 'H2O', 'CO']
        
        composition = {}
        for specie in mass_spec_list:
            composition[specie] = self.gts['%s_MS' % specie]/100.0
        ppm_list = ['C6H6','H2S','C10H8','C7H8']
        for key in ppm_list:
            composition[key] /= 10000.0
        gas_exit.pressure = self.gts.val_units('pressure_product_gas_downstream_filters')
        gas_exit.composition = composition
        gas_exit.temperature = self.gts.val_units('temp_exit_gas')
        self.gts['exit_gas_flowrate'] = gas_exit.flowrate[0]
        self.gts.units['exit_gas_flowrate'] = 'mol/s'
               

        #interpolate the Ar column to get measurements for GC comparisons
        GC_list = ['H2', 'CO', 'CO2', 'CH4', 'Ar', 'N2']
        for specie in GC_list:
            
            self.gts['%s_GC_interp' % specie] = self.gts.interpolate_col('counter', '%s_GC' % specie)

        gas_exit_GC = self.gts.outlet_stream_from_tracer([argon_tracer_feed], 'Molar', 'Ar', self.gts['Ar_GC_interp']/100.0, 'gas_exit_interp')
        gas_exit_GC.pressure = self.gts.val_units('pressure_product_gas_downstream_filters')
        gas_exit_GC.temperature = self.gts.val_units('temp_exit_gas')
        for specie in GC_list:
            gas_exit_GC.composition[specie] = self.gts["%s_GC_interp" % specie]



        self.gts.inlet_streams = [down_bed_feed, cross_brush_feed, makeup_feed, argon_tracer_feed, methane_gas_feed, steam_feed, biomass_feed]
        
        self.gts.outlet_streams = [gas_exit]
        
        self.gts.proc_elements = ['C', 'H', 'O']
        self.gts.proc_species = ['H2', 'CO', 'CO2', 'CH4', 'C2H6', 'N2', 'C2H4', 'C2H2', 'C6H6', 'Ar', 'C7H8', 'C10H8', 'C3H8', 'C3H6','H2O']
        self.gts.inert_species = ['N2', 'Ar']
        
        #set up the biomass information from the run table
        biomass_breakdown = {}
        biomass_breakdown['biomass'] = {}
        for item in ['c', 'h', 'n']:
            biomass_breakdown['biomass'][item.upper()] = self.run_info['w_%s' % item]/100.0
        biomass_breakdown['biomass']['O'] = 1 - sum(biomass_breakdown['biomass'].values())
        
        biomass_feed.special_species = biomass_breakdown
        
        

    def calculate_standard_things(self):
        """Runs the standard calculations on the data"""
        #1. Calculate the elemental and species flows
        self.gts.generate_inlet_outlet_elemental_flows()
        self.gts.generate_inlet_outlet_species_flows()

        #2. Calculate carbon conversions
        self.gts.generate_carbon_conversions()
        

        #3. Calculate changes in enthalpy and entropy
        #self.gts.generate_enthalpy_change(self, 'kW')
        #self.gts.generate_entropy_change(self, 'kW/K')

        #4. Normalize the compositions
        self.gts.generate_normalized_compositions()
        
        #calculate inlet partial pressures
        self.gts['pp_CO2'] = self.gts['CO2_inlet']/(self.gts['CO2_inlet']+self.gts['Ar_inlet']+self.gts['H2O_inlet']+self.gts['Ar_inlet'])*self.gts['pressure_reactor_gas_inlet']
        self.gts['pp_H2O'] = self.gts['H2O_inlet']/(self.gts['CO2_inlet']+self.gts['Ar_inlet']+self.gts['H2O_inlet']+self.gts['Ar_inlet'])*self.gts['pressure_reactor_gas_inlet']
        self.gts['pp_Ar'] = self.gts['Ar_inlet']/(self.gts['Ar_inlet']+self.gts['Ar_inlet']+self.gts['H2O_inlet']+self.gts['Ar_inlet'])*self.gts['pressure_reactor_gas_inlet']
        
        #5. Calculate tar loads
        self.gts.calc_tar_rate(self.gts.outlet_streams[0])

        #6. Calculate the space time
        self.gts.calc_space_time(self.reactor_size, 'biomass')

        #7. Calculate integral measures
        self.gts.generate_averages_and_stdevs(cols = self.get_analysis_config_list())
        
        #8. Calculate uncertainties
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
        
        self.gts.SQL_db_upload(self.interface_proc, table = "gas_proc_data_tbl")

        print "completed upload of timeseries data"
        #upload the integral data
        self.upload_integral_data()


        
        
    def upload_integral_data(self):
        #need to build the query elements
        objects = {}
        objects['run_id'] = str(self.run_id)
        objects['moisture'] = str(self.run_info.info['moisture'])
        objects['d10'] = str(self.run_info.info['d10'])
        objects['d50'] = str(self.run_info.info['d50'])
        objects['d90'] = str(self.run_info.info['d90'])
        objects['campaign_id'] = str(self.run_info.info['campaign_id'])
        objects['w_c'] = str(self.run_info.info['w_c'])
        objects['N_total'] = str(self.gts.numrows())
        objects['N_MS'] = str(len(self.gts['CO_MS'][np.isfinite(self.gts['CO_MS'].astype(float))]))
        for key in self.gts.avgs:
            objects['%s_avg' % key] = str(self.gts.avgs[key])
        for key in self.gts.stdevs:
            objects['%s_std' % key] = str(self.gts.stdevs[key])
        
        #try to upload the data
        try:
            q = db.insert_Query(objects, table = "gas_integral_tbl")
            #print q.getQuery()
            self.interface_proc.query(q)
        except db.DBToolboxError:
            try:
                
                q = db.update_Query(objects, table = "gas_integral_tbl", conditions = ["run_id = %s" % objects['run_id']])
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


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = "Run a gasifier analysis")
    parser.add_argument('--run_range', type=str, action = 'store')
    parser.add_argument('--run_id',type=int, action ='store')
    args = parser.parse_args()

    print args.run_range
    print args.run_id

    if args.run_id is not None:
        run_id_list = [args.run_id]

    #should really use a mutually exclusive group here, but no biggie
    elif args.run_range is not None:
        #need to parse the list -- should do some usage checking here
        main_list = args.run_range.split(",")
        run_id_list = []
        for sublist in main_list:
            if ":" in sublist:
                left = sublist.split(":")[0]
                right = sublist.split(":")[1]
                run_id_list.extend(range(int(left), int(right)+1))
            else:
                run_id_list.append(int(sublist))
        

    for run_id in run_id_list:
        print "Analyzing run %s..." % run_id
        analyzer = GasifierDataAnalysis(run_id = run_id)
        print "Data loaded"
        analyzer.calculate_standard_things()
        print "Standard things calculated"
        #analyzer.generate_output_file('output1.csv')
        analyzer.upload_to_database()
        print "Data uploaded to database"
