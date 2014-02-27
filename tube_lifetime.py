import MySQLdb as db
import argparse
import getpass
import pandas as pd
import datetime


class TubeLifetimeReport:
    """The basic data analysis class for gasifier experiments"""

    def __init__(self, user, password):
        
        tube_removal_dates = (
                    datetime.datetime(*(2012,1,1)),
                    datetime.datetime(*(2012,4,19)),
                    datetime.datetime(*(2012,5,25)),
                    datetime.datetime(*(2012,7,19)),
                    datetime.datetime(*(2012,8,22)),
                    datetime.datetime(*(2012,9,18)),
                    datetime.datetime(*(2012,11,13)),
                    datetime.datetime(*(2012,11,14)),
                    datetime.datetime(*(2013,1,17)),
                    datetime.datetime(*(2013,1,25)),
                    datetime.datetime(*(2013,2,14)),
                    datetime.datetime(*(2013,3,21)),
                    datetime.datetime(*(2013,3,27)),
                    datetime.datetime(*(2013,5,30)),
                    datetime.datetime(*(2013,6,3)),
                    datetime.datetime(*(2013,8,2)),
                    datetime.datetime(*(2014,2,4)),
                    datetime.datetime(*(2014,2,20)),
                    datetime.datetime.today()
                    )
                           
        conn = db.connect(user = user, host = '192.168.13.51', passwd = password, db = "lab_run_db")
        report = pd.DataFrame(columns = ['Start', 'End', 'Minutes Hot', 'Minutes Hot with Biomass', 'Minutes Hot with Steam', 'Minutes Hot with Biomass and Steam'])
        for i in range(len(tube_removal_dates)-1):
            start = datetime.datetime.strftime(tube_removal_dates[i], "%Y-%m-%d %H:%M:%S")
            end = datetime.datetime.strftime(tube_removal_dates[i+1], "%Y-%m-%d %H:%M:%S")
            print "This could take a minute..."
            print "Investigating tube life from %s to %s\n" %(start, end)
            df = pd.io.sql.read_frame("SELECT ts, TI_510, fr_q, MF_201, FT_310_SP FROM lv_data_gas_tbl \
                                                WHERE ts BETWEEN '%s' AND '%s' \
                                                AND TI_510>800 AND TI_510<1600 \
                                                AND ts LIKE '%%:00'" %(start, end), conn)
            
            d = pd.Series({'Start':start,
                 'End':end,
                 'Minutes Hot':len(df.index), 
                 'Minutes Hot with Biomass':len(df.query('(fr_q>0.3 & fr_q<10) | (MF_201>0.3 & MF_201<10)').index), 
                 'Minutes Hot with Steam':len(df.query('(FT_310_SP>1)').index), 
                 'Minutes Hot with Biomass and Steam':len(df.query('((fr_q>0.3 & fr_q<10) | (MF_201>0.3 & MF_201<10)) & (FT_310_SP>1)').index)})
            report.loc[i+1] = d    
            
        report.to_csv('tube_lifetime.csv')
        conn.close()
if __name__ == '__main__':

    user = raw_input('User: ')
    pswd = getpass.getpass()
    
    print "Beginning tube lifetime log.\n"
    
    tlr = TubeLifetimeReport(user, pswd)
    
