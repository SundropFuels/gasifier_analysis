from LabFunctionLib import *
from plots_toolbox import *
import dataFrame_v2
import datetime
import numpy as np
import db_toolbox as db
import argparse

class ssFinder:

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
        self._load_run_info()
        self._load_ts_timeseries_data()
        
        plot = TimeSeriesPlot(data = self.ts, Y_cols = [['mass_flow_brush_feeder', 'CO_MS', 'H2_MS']], markers='o')
        plot.plot()
        plot.show()
        plot.close()
        
    def _load_run_info(self):
        self.run_info.SQL_load(self.interface_proc, table = 'gas_run_info_tbl', run_id = self.run_id)
        
    def _load_ts_timeseries_data(self):
        """Loads raw gasifier data."""
        self.ts = GasifierProcTS(start = self.run_info.info['ts_start'], end = self.run_info.info['ts_stop'])
        self.ts.SQL_load(self.interface_raw,'gasifier_lv_GC_view') #This line needs to automatically load the units
        
        #Need to build the glossary using the SQL tools
        q = db.select_Query(objects = ['tag_number', 'simple_name', 'units'], table = "tag_glossary_tbl")
        glossary = self.interface_raw.query(q)
        self.glossary = {}
        self.gl_units = {}
        for row in glossary:
            self.glossary[row['tag_number']] = row['simple_name']
            self.gl_units[row['simple_name']] = row['units']
        
        self.ts.glossary_replace(self.glossary)
        self.ts.set_units(self.gl_units)
        self.ts.replace_None_w_nan_all()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = "Find Steady State Estimates")
    parser.add_argument('--run_range', type=str, action = 'store')
    parser.add_argument('--run_id',type=int, action ='store')
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
        print "Generating a Plot for Run %s..." % run_id
        ssfinder = ssFinder(run_id = run_id)
