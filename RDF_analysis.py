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
        self.interface_raw = db.db_interface(host = "192.168.10.20", user = "chris", passwd = "cmp87ud01")
        self.interface_raw.connect()
        q = db.use_Query("sunulator2")
        self.interface_raw.query(q)

        self.interface_proc = db.db_interface(host = "192.168.10.20", user = "chris", passwd = "cmp87ud01")
        self.interface_proc.connect()
        q = db.use_Query("sunulator2")
        self.interface_proc.query(q)


        self.run_id = run_id

        

        self.run_info = RunInformation()
        
        self._load_run_info()
        
        self._load_timeseries_data()
        self._setup_standard_streams()
        

    def _load_run_info(self):

        
        #see if we can get a complete view from this
        self.run_info.SQL_load(self.interface_proc, table = 'run_summary', run_id = self.run_id)
        #set up the tube size
        self.reactor_size = (np.pi/4*3.5*3.5*2*39.37, 'in^3')
         
               
    def _load_timeseries_data(self):
        """Loads the timeseries data into the database"""
        
        self.gts = GasifierProcTS(start = self.run_info.info['ss_start'], end = self.run_info.info['ss_end'])
        self.gts.SQL_load(self.interface_raw,'merged_data', glossary = 'glossary') #This line needs to automatically load the units
        #Need to build the glossary using the SQL tools
        q = db.select_Query(objects = ['tag_number', 'simple_name', 'units'], table = "glossary")
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
        #Need five streams for the RDF
        # 1) Biomass feed
        # 2) Entrainment gas feed (which we will adjust for the degasser)
        # 3) Sweep gas N2 feed
        # 4) Argon feed
        # 5) Steam feed
        
        #Adjust the biomass feed for RPM
        self.gts['biomass_feedrate'] = self.gts['massflowrate'] * self.run_info['feeder_slope'] + self.run_info['feeder_intercept']
        self.gts.units['biomass_feedrate'] = self.gts.units['massflowrate']

        biomass_feed = Stream('biomass_feed',flowrate = self.gts.val_units('biomass_feedrate'),composition = {'H2O':self.run_info['biomass_moisture']/100.0, 'biomass':1.00-self.run_info['biomass_moisture']/100.0}, basis = "mass")

        #create a new entrainment gas column
        
        self.gts['entrainment'] = self.gts['degasser_active'] * self.gts['net_entrainment_flow'] + (np.ones(len(self.gts['degasser_active']))-self.gts['degasser_active'])*self.gts['ent_gas_flow']
        if (self.gts['degasser_active'] == 1).any():
            self.gts.units['entrainment'] = 'L/min'
        else:
            self.gts.units['entrainment'] = 'ft^3/min'
        
        entrainment_gas = Stream('entrainment_feed', flowrate = self.gts.val_units('entrainment'), composition = {'N2':1.0}, basis = "std_gas_volume") 

        sweep_gas = Stream('sweep_feed', flowrate = self.gts.val_units('sweep_flow'), composition = {'N2':1.0}, basis = "std_gas_volume")
        
        tracer_gas = Stream('tracer_feed', flowrate = self.gts.val_units('tracer_flow'), composition = {'Ar':1.0}, basis = "std_gas_volume")

        steam_feed = Stream('steam_feed', flowrate = self.gts.val_units('steamflowrate'), composition = {'H2O':1.0}, basis = "mass")       
                     
        #Set stream pressures and temperatures
        MFC_SP = (101325.0, 'Pa')
        MFC_ST = (70.0, 'F')

        biomass_feed.set_temperature((25.0, 'C'))
        entrainment_gas.set_temperature((25.0, 'C'))
        sweep_gas.set_temperature((25.0, 'C'))
        tracer_gas.set_temperature((25.0, 'C'))
        steam_feed.set_temperature(self.gts.val_units('steamtemp'))
        
        biomass_feed.set_pressure(self.gts.val_units('pressure'))
        entrainment_gas.set_pressure(self.gts.val_units('pressure'))
        sweep_gas.set_pressure(self.gts.val_units('pressure'))
        tracer_gas.set_pressure(self.gts.val_units('pressure'))
        steam_feed.set_pressure(self.gts.val_units('pressure'))
               
        
        #Setup the standard temperatures and pressures
        entrainment_gas.std_temperature = MFC_ST
        sweep_gas.std_temperature = MFC_ST
        tracer_gas.std_temperature = MFC_ST
        
        entrainment_gas.std_pressure = MFC_SP
        sweep_gas.std_pressure = MFC_SP
        tracer_gas.std_pressure = MFC_SP
        
        #Set up exit gas flowrates
        
        
        gas_exit = self.gts.outlet_stream_from_tracer([tracer_gas],"Molar", "Ar", self.gts['mol_pct_Ar']/100.0, 'gas_exit')  #may need to consider using interpolated data -- how does this work out for the current experiments?
        #Need to set up the exit gas composition now -- we can actually build a variety of different streams here and use them as necessary 

        gc_list = ['H2', 'O2', 'N2', 'CH4', 'CO', 'CO2', 'Ar', 'C2H6', 'C2H4', 'C2H2', 'C3H8', 'C3H6']
                
        composition = {}
        for specie in gc_list:
            composition[specie] = self.gts['mol_pct_%s' % specie]/100.0
        
        gas_exit.set_pressure(self.gts.val_units('pressure'))
        gas_exit.set_composition(composition)
        gas_exit.set_temperature(self.gts.val_units('producttemp'))
        self.gts['exit_gas_flowrate'] = gas_exit.flowrate[0]
        self.gts.units['exit_gas_flowrate'] = 'mol/s'
               

        #create an interpolated column -- may be necessary if we have short runs/very sparse data
        
        for specie in gc_list:
            
            self.gts['mol_pct_%s_interp' % specie] = self.gts.interpolate_col('counter', 'mol_pct_%s' % specie)

        gas_exit_interp = self.gts.outlet_stream_from_tracer([tracer_gas], 'Molar', 'Ar', self.gts['mol_pct_Ar_interp']/100.0, 'gas_exit_interp')
        gas_exit_interp.set_pressure(self.gts.val_units('pressure'))
        gas_exit_interp.set_temperature(self.gts.val_units('producttemp'))
        gc_composition = {}
        for specie in gc_list:
            gc_composition[specie] = self.gts["mol_pct_%s_interp" % specie]
            gas_exit_interp.set_composition(gc_composition)


        self.gts.inlet_streams = [biomass_feed, entrainment_gas, tracer_gas, sweep_gas, steam_feed]
        
        self.gts.outlet_streams = [gas_exit]
        
        self.gts.proc_elements = ['C', 'H', 'O']
        self.gts.proc_species = ['H2', 'CO', 'CO2', 'CH4', 'C2H6', 'N2', 'C2H4', 'C2H2', 'Ar', 'C3H8', 'C3H6','H2O']
        self.gts.inert_species = ['N2', 'Ar']
        
        #set up the biomass information from the run table
        biomass_breakdown = {}
        biomass_breakdown['biomass'] = {}
        for item in ['C', 'H', 'N']:
            biomass_breakdown['biomass'][item] = self.run_info['%s_pct' % item]/100.0
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
        self.gts.generate_enthalpy_change('kW')
        #self.gts.generate_entropy_change(self, 'kW/K')

        #4. Normalize the compositions
        self.gts.generate_normalized_compositions()
        
        #calculate inlet partial pressures
        self.gts['pp_CO2'] = self.gts['CO2_inlet']/(self.gts['CO2_inlet']+self.gts['Ar_inlet']+self.gts['H2O_inlet']+self.gts['Ar_inlet']+self.gts['N2_inlet'])*self.gts['pressure']  
        self.gts['pp_H2O'] = self.gts['H2O_inlet']/(self.gts['CO2_inlet']+self.gts['Ar_inlet']+self.gts['H2O_inlet']+self.gts['Ar_inlet']+self.gts['N2_inlet'])*self.gts['pressure']
        self.gts['pp_Ar'] = self.gts['Ar_inlet']/(self.gts['Ar_inlet']+self.gts['Ar_inlet']+self.gts['H2O_inlet']+self.gts['Ar_inlet']+self.gts['N2_inlet'])*self.gts['pressure']
        
        #5. Calculate tar loads - NA FOR RDF DATA
        #self.gts.calc_tar_rate(self.gts.outlet_streams[0])

        #6. Calculate the space time
        print "calculating space time..."
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
        f.write('End time: %s\n' % self.run_info['ss_end']) #May need to do some formatting here
        #f.write('Biomass type', self.run_info.biomass_type)
        #f.write('Biomass moisture', self.run_info.moisture)

        #NEED TO INCLUDE SETPOINT INFORMATION HERE, POSSIBLY

        #now, write the data to the csv file
        f.close()
        #self.gts.write_csv(filename, mode = 'append', col_list = ["timestamp", "X_tot", "X_good", "X_std", "CO2_normalized", "CO_normalized", "CH4_normalized", "H2_normalized", "C2H2_outlet", "C2H4_outlet", "CH4_outlet", "C6H6_outlet", "tar_loading", "tar_loading_incl"])
        self.gts.write_csv(filename, mode = 'append')

    def upload_to_database(self):
        
        

        
        #upload the time series data
        
        self.gts.SQL_db_upload(self.interface_proc, table = "RDF_proc_data")

        print "completed upload of timeseries data"
        #upload the integral data
        self.upload_integral_data()


        
        
    def upload_integral_data(self):
        #need to build the query elements
        objects = {}
        objects['run_id'] = str(self.run_id)
        
        for key in self.gts.avgs:
            objects['%s_avg' % key] = str(self.gts.avgs[key])
        for key in self.gts.stdevs:
            objects['%s_std' % key] = str(self.gts.stdevs[key])
        
        #try to upload the data
        try:
            q = db.insert_Query(objects, table = "RDF_integral_tbl")
            #print q.getQuery()
            self.interface_proc.query(q)
        except db.DBToolboxError:
            try:
                
                q = db.update_Query(objects, table = "RDF_integral_tbl", conditions = ["run_id = %s" % objects['run_id']])
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
    parser.add_argument('--run_id',type=str, action ='store')
    parser.add_argument('--file',type=str,action = 'store')
    args = parser.parse_args()

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
        analyzer = GasifierDataAnalysis(run_id = run_id)
        #print "Data loaded"
        analyzer.calculate_standard_things()
        #print "Standard things calculated"
        #analyzer.generate_output_file('RDF007-027.csv')
        analyzer.upload_to_database()
        #print "Data uploaded to database"



