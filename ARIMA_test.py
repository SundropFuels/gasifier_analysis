#!/usr/bin/env python

from LabFunctionLib import *
import datetime
import numpy as np
import db_toolbox as db
import Thermo
import argparse
import csv
import dataFrame_v2 as df
import statsmodels.tsa.arima_model as ARIMA


class ARIMA_test:
    """The basic data analysis class for gasifier experiments"""

    def __init__(self, run_id, user, password, run_information = None):
        #Create the gasifier data frame, and load the data from the SQL database (this will be automated through the interface later)
        
        self.interface_proc = db.db_interface(host = "192.168.13.51", user = user, passwd = password)
        self.interface_proc.connect()
        q = db.use_Query("lab_proc_db")
        self.interface_proc.query(q)
        self.run_id = run_id
        
        self.run_info = RunInformation()
        #self.run_info.info = run_information
        self._load_run_info()
        self._load_timeseries_data()
        

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


    def _load_timeseries_data(self):

        self.gts = df.Dataframe()
        self.gts.SQL_load_data(self.interface_proc,'gas_proc_data_tbl', conditions = ["timestamp >= '%s'" % self.run_info.info['ss_start'],"timestamp < '%s'" % self.run_info.info['ss_stop']]) #This line needs to automatically load the units
        

if __name__ == "__main__":

    user = raw_input('User: ')
    pswd = getpass.getpass()
    
    test = ARIMA_test(run_id = 106, user = user, password = pswd)
    
    print test.gts['timestamp']
    model = ARIMA.ARMA(test.gts['mass_flow_brush_feeder'], order = (1,1))
    result = model.fit()
    print result.summary()
    print test.gts['mass_flow_brush_feeder']
    print result.fittedvalues
    print result.resid
    
	
