#!/usr/bin/env python

import LabFunctionLib as lfl
import datetime
import numpy as np
import db_toolbox as db
import argparse
import getpass
import pandas as pd

class logException(Exception):
    def __init__(self,value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class timestampError(logException):
    pass

class log:
    def __init__(self, start, stop):
       
        if type(start) ==  datetime.datetime: 
            self.start = start
        else: 
            try:
                self.start = datetime.datetime.strptime(start, "%Y-%m-%d %H:%M:%S")
            except:
                raise timestampError, "Start time must either be datetime.datetime object or in 'YYYY-mm-dd H:M:S' format"
                
        if type(stop) ==  datetime.datetime:
            self.stop = stop
        else: 
            try:
                self.stop = datetime.datetime.strptime(stop, "%Y-%m-%d %H:%M:%S")
            except:
                raise timestampError, "Stop time must either be datetime.datetime object or in 'YYYY-mm-dd H:M:S' format"
        
    def SQL_connect(self, host, user, password, database):
        """Creates database interface for connecting to SQL database"""
        
        self.db_interface = db.db_interface(host = host, user = user, passwd = password)
        self.db_interface.connect()
        q = db.use_Query(database)
        self.db_interface.query(q)
             
    def load_SQL_time_series(self, con = None, table = ""):
        """Loads time series data for specified time frame"""
        
        if con == None:
            con = self.db_interface
        self.ts = lfl.ts_data(self.start, self.stop)
        q = db.select_Query(objects = ['tag_number'], table = "tag_glossary_tbl", condition_list = ["tag_number LIKE 'MFC_%_SP'"])
        self.log_tags = [i['tag_number'] for i in con.query(q)]
        self.ts.SQL_load(con, table = table, glossary = 'tag_glossary_tbl')
        self.log_status = {}
        self.glossary = {}
        self.gl_units = {}
        
        for row in self.glossary:
            self.glossary[row['tag_number']] = row['simple_name']
            self.gl_units[row['simple_name']] = row['units']
            self.log_status[row['log_status']] = row['log_status']
        
        self.ts.glossary_replace(self.glossary)
        self.ts.set_units(self.gl_units)
        self.ts.replace_None_w_nan_all()

    def find_changes(self, column_name):
        """Builds a list of index positions for a column where the values change"""
        col = self.ts[column_name]
        changes = []
        
        i= 0
        while i < len(col):      
            try: 
                s = col[i:]
                i = s[s!=s[i]].index[0]
                if i == 0: break
                changes.append(i)              
            except IndexError:
                break
           
        return changes       
        
    def generate_log_dataframe_(self):
        """Generates Pandas DataFrame for system set point changes from time series data"""
        self.log_df = pd.DataFrame(columns = ['ts', 'tag_number', 'old_value', 'new_value'])
        for col in self.log_tags:
            changes = self.find_changes(col)
            for i in changes:
                self.log_df = self.log_df.append({'ts':self.ts['ts'][i], 
                                                  'tag_number':col, 
                                                  'old_value':self.ts[col][i-1], 
                                                  'new_value':self.ts[col][i]},
                                                 ignore_index = True)
        self.log_df = self.log_df.sort(column = 'ts')                                         
        print self.log_df
    def upload_log_dataframe(self, con = None, table = 'system_log_tbl'):
        """Uploads log dataframe onto database."""
        
        if con == None:
            con = self.db_interface
        
        for i in self.log_df.index:
            try:
                self.log_df.loc[i].to_sql(table, con, flavor = 'mysql', if_exitst = 'append')
            except:
                pass
        
    def write_log_book(self, con = None, table = 'system_log_tbl'):
        """Writes human readable log book from system changes and operator comments"""
        pass
        

        
        
        
