import dataFrame_v2 as df
import datetime
import numpy as np
import matplotlib.pyplot as plt
import db_toolbox as db
import argparse

class TubeHistoryReport:

    def __init__(self, start, stop):
        #Create the gasifier data frame, and load the data from the SQL database
        self.start=start
        self.stop=stop
        self.interface_proc = db.db_interface(host = "192.168.10.20", user = "chris", passwd = "cmp87ud01")
        self.interface_proc.connect()
        q = db.use_Query("lab_proc_db")
        self.interface_proc.query(q)
        self._load_dataframe()
        self.tubehistorysummary()
        self.write_csv()

    def _load_dataframe(self):
        self.data=df.Dataframe()
        self.data.SQL_load_data(self.interface_proc, "gas_run_info_tbl", ["ts_start >= '%s' and ts_start < '%s'" %(self.start,self.stop)])
        self.data['hours']=np.array([i.seconds/3600. for i in (self.data['ts_stop']-self.data['ts_start'])])

    def tubehistorysummary(self):
        temperatures=np.unique(self.data['tic510_sp'])
        feedrates=np.unique(self.data['biomass_rate'])
        self.summary=[]
        for temp in temperatures:
            for fr in feedrates:
                time=0
                for run in range(len(self.data['run_id'])):
                    if self.data['tic510_sp'][run]==temp and self.data['biomass_rate'][run]==fr:
                        time+=self.data['hours'][run]
                self.summary.append({'temperature':temp, 'biomass feedrate':fr, 'hours':time})
        print 'Tube run history for %s through %s' %(self.start, self.stop)
        print 'Temperature\tBiomass Feedrate\tHours'
        for i in range(len(self.summary)):
            if self.summary[i]['hours']>0:
                print self.summary[i]['temperature'], '\t\t', self.summary[i]['biomass feedrate'], '\t\t\t', self.summary[i]['hours']
                
    def write_csv(self):
        with open('tube_history.csv','a') as f:
            f.write('Tube run history for %s through %s\n' %(self.start, self.stop))
            f.write('Temperature,Biomass Feedrate,Hours\n')
            for i in range(len(self.summary)):
                if self.summary[i]['hours']>0:
                    f.write(str(self.summary[i]['temperature'])+',')
                    f.write(str(self.summary[i]['biomass feedrate'])+',')
                    f.write(str(self.summary[i]['hours'])+'\n')
            f.write('\n')
                
        
                                    
                    
if __name__=='__main__':
    parser = argparse.ArgumentParser(description = "Run a gasifier analysis")
    parser.add_argument('--start', type=str, action = 'store')
    parser.add_argument('--stop',type=str, action ='store')
    args = parser.parse_args()

    if args.start is not None:
        start = args.start
        
    if args.stop is not None:
        stop = args.stop

    #start='2012-12-10 00:00'
    #stop='2013-02-28 00:00'
    tubehistory=TubeHistoryReport(start,stop)
    
