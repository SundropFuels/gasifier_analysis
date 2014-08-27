import MySQLdb as db
import argparse
import getpass
import pandas as pd
import datetime
import numpy as np

class TubeLifetimeReport:
    """The basic data analysis class for gasifier experiments"""

    def __init__(self, user, password):
        
        tube_removal_dates = (
                    datetime.datetime(*(2012,1,1)),
                    datetime.datetime(*(2012,4,19,17,0)),
                    datetime.datetime(*(2012,5,25)),
                    datetime.datetime(*(2012,7,19,18,0)),
                    datetime.datetime(*(2012,8,22)),
                    datetime.datetime(*(2012,9,18)),
                    datetime.datetime(*(2012,11,13,5,0)),
                    datetime.datetime(*(2012,11,15,20,0)),
                    datetime.datetime(*(2013,1,16,8,0)),
                    datetime.datetime(*(2013,1,17,1,0)),
                    datetime.datetime(*(2013,1,25,21,0)),
                    datetime.datetime(*(2013,2,14,6,0)),
                    datetime.datetime(*(2013,3,21)),
                    datetime.datetime(*(2013,3,27,6,0)),
                    datetime.datetime(*(2013,5,30,15,0)),
                    datetime.datetime(*(2013,6,3)),
                    datetime.datetime(*(2013,8,2,11,0)),
                    datetime.datetime(*(2014,2,4,18,0)),
                    datetime.datetime(*(2014,2,20,13,0)),
                    datetime.datetime.today()
                    )
                           
        conn = db.connect(user = user, host = '192.168.13.51', passwd = password, db = "lab_run_db")
        report = pd.DataFrame(columns = ['Start', 'End', 'Minutes Hot', 'Minutes Hot with Biomass', 'Minutes Hot with Steam', 'Minutes Hot with Biomass and Steam',
                                         'Total Biomass Fed (lbs)', 'Total Steam Fed (lbs)', 'Minutes 800-1000 C', 'Minutes 1000-1200 C', 'Minutes 1200-1300 C', 'Minutes 1300-1400 C',
                                         'Minutes 1400-1500 C', 'Minutes Hot 0-20 psig', 'Minutes Hot 20-40 psig', 'Minutes Hot > 40 psig'])
        for i in range(len(tube_removal_dates)-1):
            start = datetime.datetime.strftime(tube_removal_dates[i], "%Y-%m-%d %H:%M:%S")
            end = datetime.datetime.strftime(tube_removal_dates[i+1], "%Y-%m-%d %H:%M:%S")
            print "This could take a minute..."
            print "Investigating tube life from %s to %s\n" %(start, end)
            
            # Query looks at one point at 0 seconds every minute and sees if it is at temperature.  Length of dataframes is assumed to be length at temperature in minutes.
            df = pd.io.sql.read_frame("SELECT ts, TI_510, fr_q, MF_201, FT_310_SP, PI_600 FROM lv_data_gas_tbl \
                                                WHERE ts BETWEEN '%s' AND '%s' \
                                                AND TI_510>800 AND TI_510<1600 \
                                                AND ts LIKE '%%:00'" %(start, end), conn)
            
            len_brush=0
            avg_brush=1
            len_ktron=0
            avg_ktron=1
            
            if len(df.query('(fr_q>0.3 & fr_q<10)')) >0:
                len_brush = len(df.query('(fr_q>0.3 & fr_q<10)').index)
                avg_brush = df.query('(fr_q>0.3 & fr_q<10)')['fr_q'].mean()
            if len(df.query('(MF_201>0.3 & MF_201<10)')) >0:
                len_brush = len(df.query('(MF_201>0.3 & MF_201<10)').index)
                avg_brush = df.query('(MF_201>0.3 & MF_201<10)')['MF_201'].mean()    
            
            
            d = pd.Series({'Start':start,
                 'End':end,
                 'Minutes Hot':len(df.index), 
                 'Minutes Hot with Biomass':len(df.query('(fr_q>0.3 & fr_q<10) | (MF_201>0.3 & MF_201<10)').index), 
                 'Minutes Hot with Steam':len(df.query('(FT_310_SP>1)').index), 
                 'Minutes Hot with Biomass and Steam':len(df.query('((fr_q>0.3 & fr_q<10) | (MF_201>0.3 & MF_201<10)) & (FT_310_SP>1)').index),
                 'Total Biomass Fed (lbs)': round(avg_brush/60.*len_brush + avg_ktron/60.*len_ktron,2), #Average biomass feed rate of column times length of column
                 'Total Steam Fed (lbs)': round(df['FT_310_SP'].mean()*len(df.index)*0.00220462,2), #Average steam feed rate of column times length of column 
                 'Minutes 800-1000 C':len(df.query('((TI_510 > 800) & (TI_510 < 1000))').index),
                 'Minutes 1000-1200 C':len(df.query('((TI_510 > 1000) & (TI_510 < 1200))').index),
                 'Minutes 1200-1300 C':len(df.query('((TI_510 > 1200) & (TI_510 < 1300))').index),
                 'Minutes 1300-1400 C':len(df.query('((TI_510 > 1300) & (TI_510 < 1400))').index),
                 'Minutes 1400-1500 C':len(df.query('((TI_510 > 1400) & (TI_510 < 1500))').index),
                 'Minutes Hot 0-20 psig':len(df.query('((PI_600 > 0) & (PI_600 < 20))').index), 
                 'Minutes Hot 20-40 psig':len(df.query('((PI_600 > 20) & (PI_600 < 40))').index), 
                 'Minutes Hot > 40 psig':len(df.query('((PI_600 > 40))').index)
                 })
            report.loc[i+1] = d    
            
        report.to_csv('tube_lifetime.csv')
        conn.close()
if __name__ == '__main__':

    user = raw_input('User: ')
    pswd = getpass.getpass()
    
    print "Beginning tube lifetime log.\n"
    
    tlr = TubeLifetimeReport(user, pswd)
    
