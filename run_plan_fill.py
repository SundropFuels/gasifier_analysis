"""This is the CSV uploading code for run sheets.  It allows one to use e.g. Excel, Calc to act as a GUI"""

import pandas as pd
import db_toolbox as db
import numpy as np
import argparse
import getpass
import pandas.io.sql as psql

class PropertyRow:
    """This is the abstract class for various conditions (setpoints, biomass, etc)"""
    def __init__(self, **kwargs):

        self.properties = {}

        for key in kwargs:
            self.properties[key] = kwargs[key]


    def __eq__(self, other):
        if not isinstance(other, type(self)):
            raise TypeError, "Cannot compare type %s with PropertyTable" % type(other)

        flag = True
        for key in self.properties:
            if self.properties[key] != other.properties[key]:
                flag = False

        return flag

    def __ne__(self, other):
        return not __eq__(self, other)


class SetpointRow(PropertyRow):

    def __init__(self, **kwargs):
        keylist = ['temperature', 'pressure', 'biomass_rate', 'steam_flow', 'steam_temp', 'ent_CO2', 'sweep_CO2', 'Ar_tracer', 'superheater_purge', 'tube_diameter']
        for key in kwargs:
            if key not in keylist:
                raise Exception, "%s is not an appropriate key for the setpoints" % key
        for key in keylist:
            if key not in kwargs:
                raise Exception, "%s is missing from the arguments" % key

        PropertyRow.__init__(self, **kwargs)

class BiomassRow(PropertyRow):

    def __init__(self, **kwargs):
        keylist = ['sample_name', 'moisture', 'w_c', 'w_n', 'w_h', 'd10', 'd50', 'd90']
        for key in kwargs:
            if key not in keylist:
                raise Exception, "%s is not an appropriate key for the setpoints" % key
        for key in keylist:
            if key not in kwargs:
                raise Exception, "%s is missing from the arguments" % key

        PropertyRow.__init__(self, **kwargs)
    

class RunTableUploader:

    def __init__(self, user, password, host = 'localhost'):
        
        #set up the database connection
        self.interface = db.db_interface(host, user = user, passwd = password)
        self.interface.connect()
        q = db.use_Query("pilot_proc_db")
        self.interface.query(q)


        #create any appropriate member objects

        pass

    def read_data(self, filename):
        #load in the data
        self.data = pd.read_csv(filename)

        self.data['setpoint_id'] = np.zeros(len(self.data.index))
        self.data['biomass_id'] = np.zeros(len(self.data.index))

        #Now, need to fill in the setpoint_id and biomass_id if they exist
        #Create a unique set of setpoint_ids and biomass_ids from the pandas dataframe
        setpoints = []
        biomasses = []
        
        sp_info = ['temperature', 'pressure', 'biomass_rate', 'steam_flow', 'steam_temp', 'ent_CO2', 'sweep_CO2', 'Ar_tracer', 'superheater_purge', 'tube_diameter']

        bm_info = ['sample_name', 'moisture', 'w_c', 'w_n', 'w_h', 'd10', 'd50', 'd90']

        for i in self.data.index:
            d = {}  
            for j in sp_info:
                d[j] = self.data[j][i]
            setp = SetpointRow(**d)
            d = {}
            for j in bm_info:
                d[j] = self.data[j][i]
            bmp = BiomassRow(**d)
            spf = False
            bmf = False
            for sp in setpoints:
                if setp == sp:
                    spf = True
            for bm in biomasses:
                if bmp == bm:
                    bmf = True
            if not spf:
                setpoints.append(setp)
            if not bmf:
                biomasses.append(bmp)

        #Now we need to get the existing setpoint and biomass tables
        try:
            setpoint_table = pd.read_sql("SELECT * FROM setpoint_tbl", con = self.interface.connection)
        
        except AttributeError:
            setpoint_table = psql.read_frame("SELECT * FROM setpoint_tbl", con = self.interface.connection)


        try:
            biomass_table = pd.read_sql("SELECT * FROM biomass_tbl", con = self.interface.connection)

        except AttributeError:
	    biomass_table = psql.read_frame("SELECT * FROM biomass_tbl", con = self.interface.connection)

        #Check and see if a specific setpoint is already in the table
        for i in setpoint_table.index:
            d = {}
            for j in sp_info:
                d[j] = setpoint_table[j][i]
            setp = SetpointRow(**d)

        for i in biomass_table.index:
            d = {}
            for j in bm_info:
                d[j] = biomass_table[j][i]
            bmp = BiomassRow(**d)

            for s in setpoints:
                if s == setp:
                    setpoints.remove(s)
                    break #all the items in setpoint are already unique

            for b in biomasses:
                if b == bmp:
                    biomasses.remove(b)
                    break #all the items in biomasses are already unique

        self.setpoints = setpoints
        self.biomasses = biomasses
        
       

    def append_runsheet(self):
        #upload the unique setpoint and biomass points (which will auto-increment ids)
        for setp in self.setpoints:
            objects = {}
            for k in setp.properties:
                #print setp.properties[k]
                objects[k] = str(setp.properties[k])

            
            q = db.insert_Query(objects = objects, table = "setpoint_tbl")
            #print q.getQuery()
            self.interface.query(q)

        for bm in self.biomasses:
            objects = {}
            for k in bm.properties:
                #print bm.properties[k]
                objects[k] = str(bm.properties[k])
            
            q = db.insert_Query(objects = objects, table = "biomass_tbl")
            #print q.getQuery()
            self.interface.query(q)
        #now upload the run information, changing to NULL where appropriate

        

        run_info = ['run_id', 'exp_id', 'ts_start', 'ts_stop', 'ss_start', 'ss_stop', 'operator', 'feeder_slope', 'feeder_intercept']

        sp_info = ['temperature', 'pressure', 'biomass_rate', 'steam_flow', 'steam_temp', 'ent_CO2', 'sweep_CO2', 'Ar_tracer', 'superheater_purge', 'tube_diameter']

        bm_info = ['sample_name', 'moisture', 'w_c', 'w_n', 'w_h', 'd10', 'd50', 'd90']

        for i in self.data.index:
            d = {}
            for k in run_info:
                d[k] = str(self.data[k][i])
            #now locate the setpoint in the setpoint_table using SQL queries
            objects = ['setpoint_id']
            
            cl = []
            for item in sp_info:
                cl.append("%s='%s'" % (item, self.data[item][i]))
            q = db.select_Query(objects = objects, condition_list = cl, table = "setpoint_tbl")
            s = self.interface.query(q)
            sid = s[0]['setpoint_id']

            cl = []
            objects = ['biomass_id']
            for item in bm_info:
                cl.append("%s='%s'" % (item, self.data[item][i]))
            q = db.select_Query(objects = objects, condition_list = cl, table = "biomass_tbl")
            s = self.interface.query(q)
            bid = s[0]['biomass_id']

            d['setpoint_id'] = str(sid)
            d['biomass_id'] = str(bid)

            
            try:
                 q = db.insert_Query(objects = d, table = "run_info_tbl")
                 self.interface.query(q)
            except db.DBToolboxError:
                try:
                    q = db.update_Query(objects = d, table = "run_info_tbl", conditions = ["%s = '%s'" % ('run_id', d['run_id'])])
                    self.interface.query(q)
                except db.DBToolboxError:
                    raise Exception, "Trouble loading the data, man.  How bout you go catch some tasty waves?"         
        

    def upload_runsheet(self):
        pass


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = "Run a gasifier analysis")
    parser.add_argument('--filename', type=str, action = 'store')
    parser.add_argument('--host', type = str, action = 'store')
    parser.add_argument('--user', type = str, action = 'store')
    parser.add_argument('--pswd', type = str, action = 'store')
    args = parser.parse_args()

    if args.user is not None and args.pswd is not None:
        user = args.user
        pswd = args.pswd
    else:
        user = raw_input('User: ')
        pswd = getpass.getpass()
        
    if args.filename is not None:
        f = args.filename

    if args.host is not None:
        h = args.host
    else:
        h = 'localhost'

        
    print "Uploading run table" 
    loader = RunTableUploader(user = user, password = pswd, host = h)
    loader.read_data(f)
    loader.append_runsheet()
    print "done"
        
