from math import *
import os
from LabFunctionLib import *
import datetime
import numpy as np
import matplotlib.pyplot as plt
import db_toolbox as db
import Thermo
import argparse
import csv
import subprocess
import statsmodels.tsa.arima_model as ARIMA
import dataFrame_v2 as df
from plots_toolbox import *


class GasifierReport:
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
        self.run_integrals = RunInformation()
        #self.run_info.info = run_information
        self._load_run_info()
        self._create_file_structure()
        self._load_run_integrals()
        self._load_ts_timeseries_data()
        self._load_ss_timeseries_data()
        self._add_units_to_run_info()
        self._convert_gas_units_run_info()
        self._convert_conversions_to_percent()
        self._add_std_percent()
        self._convert_steam_flow_to_ml_min()
        
        ###TESTS OF NEW PLOT OBJECTS HERE!"""
        
        

        z = FourPlot(data = self.ss, x_label = 'Timestamp', y_label = 'Mass Feedrate (lb/hr)', x_var = 'timestamp', y_var = 'mass_flow_brush_feeder')
        z.plot()
        z.show()
        z.close()


        h = ControlChartfromDataframe(data = self.ss, y_col = 'mass_flow_brush_feeder', x_col = 'timestamp', sample_size = 1)
        a = h.getDataframe()

        c = IndividualsXBarControlChart(data = a, y_label = 'Mass Flow Brush Feeder', x_label = 'Time')
        c.plot()
        c.show()
        c.close()
        
        self.gas_comp_pie_plot()
        
        """
        self.gas_comp_pie_plot()
        
        self.time_series_plot(['mass_flow_brush_feeder'])
        self.time_series_plot(['CO_MS','CO2_MS', 'H2_MS', 'CH4_MS','Ar_MS'],
                              ylabel='Gas Composition (% vol)',
                              caption='Time series plot for gas composition.')
        self.time_series_plot(['C2H4_MS', 'C6H6_MS', 'C7H8_MS', 'C10H8_MS'],
                              ylabel='Gas Composition (ppm)',
                              caption='Time series plot for gas composition.')
        self.time_series_plot(['temp_skin_tube_top','temp_skin_tube_middle','temp_skin_tube_bottom'],
                              ylabel='Tube Skin Temperatures ($^\circ$C)',
                              caption='Time series plot for reactor tube skin temperatures.')
        
        #Create ARIMA fits as necessary
        ARIMA_list = ['mass_flow_brush_feeder', 'temp_skin_tube_middle', 'temp_steam_reactor_entry']
        for col in ARIMA_list:
            self.fit_ARIMA(col)



        self.control_plots = [('mass_flow_brush_feeder',10),('temp_skin_tube_middle',10),('temp_steam_reactor_entry',10),('pressure_ash_knockout_vessel',10),('CO_MS',3),('CO2_MS',3),('H2_MS',3),('CH4_MS',3)]
        for p in self.control_plots:
            self.control_plot(p[0],p[1])       

         

        #Need to get these from a database table so that it is configurable
        
        self.four_plot('mass_flow_brush_feeder',
                       ylabel='Biomass Flow Rate (lbs/hr)',
                       caption='Four-plot for biomass flow rate')
        self.four_plot('temp_skin_tube_middle',
                       ylabel='Reactor Skin Middle ($^\circ$C)',
                       caption='Four-plot for reactor skin temperature.')
        self.four_plot('temp_steam_reactor_entry',
                       ylabel='Steam Temperature ($^\circ$C)',
                       caption='Four-plot for temperature of steam at reactor inlet.')
        self.four_plot('pressure_ash_knockout_vessel',
                       ylabel='Ash Knockout Pressure (psig)',
                       caption='Four-plot for ash knockout pressure.')
        self.four_plot('CO_MS',
                       ylabel='Carbon Monoxide (%vol)',
                       caption='Four-plot for carbon monoxide readings from the mass spectrometer.')
        self.four_plot('H2_MS',
                       ylabel='Hydrogen (%vol)',
                       caption='Four-plot for hydrogen readings from the mass spectrometer.')
        self.four_plot('CO2_MS',
                       ylabel='Carbon Dioxide (%vol)',
                       caption='Four-plot for carbon dioxide readings from the mass spectrometer.')
        self.four_plot('CH4_MS',
                       ylabel='Methane (%vol)',
                       caption='Four-plot for methane readings from the mass spectrometer.')

        self.generate_standard_report()
        """
        
    def _create_file_structure(self):
        self.directory='rpt/%s/' % str(self.run_info.info['run_id'])
        if not os.path.exists(self.directory):
            os.makedirs(self.directory)
        
    def _load_run_info(self):
        self.run_info.SQL_load(self.interface_proc, table = 'gas_run_info_tbl', run_id = self.run_id)

    def _load_run_integrals(self):
        self.run_integrals.SQL_load(self.interface_proc, table = 'gas_integral_tbl', run_id = self.run_id)
        for i in self.run_integrals.info:
            if i not in self.run_info.info:
                self.run_info.info[i]=self.run_integrals.info[i]

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

    def _load_ss_timeseries_data(self):
        """Loads processed steady state data including calculated columns"""
        self.ss = ts_data(start = self.run_info.info['ss_start'], end = self.run_info.info['ss_stop'])
        self.ss.SQL_load(self.interface_proc,'gas_proc_data_tbl')
        
        self.ss.replace_None_w_nan_all()

    def _add_units_to_run_info(self):
        for i in self.ts.units:
            if self.ts.units[i]=='C':
                self.run_info.info[i+'_units']=r'$^\circ'+r'$C'
            elif self.ts.units[i]=='%':
                self.run_info.info[i+'_units']='\%'
            else:
                self.run_info.info[i+'_units']=self.ts.units[i]

    def _add_std_percent(self):
        l=self.run_info.info.keys()
        for i in l:
            if i.endswith('_std'):
                try: self.run_info.info[i[:-4]+'_pstd']=self.run_info.info[i]/self.run_info.info[i[:-4]+'_avg']*100
                except KeyError:
                    pass
                except ZeroDivisionError:
                    self.run_info.info[i[:-4]+'_pstd']=0

    ##CHANGE -- USE UNIT CONVERTER
    def _convert_gas_units_run_info(self):
        gaslist=[]
        for i in self.run_info.info:
            if i.endswith('_normalized_avg'):
                gaslist.append(i.replace('_avg', ''))
        for gas in gaslist:
            self.run_info.info[gas+'_avg']*=100
            self.run_info.info[gas+'_std']*=100
            self.run_info.info[gas+'_units']='\%'
            if self.run_info.info[gas+'_avg']<1 and self.run_info.info[gas+'_avg']>0:
               self.run_info.info[gas+'_avg']*=10000
               self.run_info.info[gas+'_std']*=10000
               self.run_info.info[gas+'_units']='ppm'
            if self.run_info.info[gas.replace('_normalized','_MS_units')]=='ppm':
                self.run_info.info[gas+'_units']='ppm'                                  
        normprod=self.run_info.info['CO_normalized_avg']/self.run_info.info['CO_MS_avg']
        self.run_info.info['H2S_normalized_avg']=self.run_info.info['H2S_MS_avg']*normprod
        self.run_info.info['H2S_normalized_std']=self.run_info.info['H2S_MS_std']*normprod
        self.run_info.info['H2S_normalized_units']='ppm'
    

    ##CHANGE        
    def _convert_steam_flow_to_ml_min(self):
        self.run_info.info['steam_flow_avg']*=7.55987283

    def _convert_conversions_to_percent(self):
        l=[i for i in self.run_info.info if i.startswith('X_')]
        for i in l:
            self.run_info.info[i]*=100
        
    
      
    def gas_comp_pie_plot(self):
        gasdict={}
        for i in self.run_info.info:
            if i.endswith('_MS_avg'):
                gasdict[i]=self.run_info.info[i]
        goodgas=['CO', 'CO2', 'CH4', 'H2']
        targas=['C2H6', 'C2H4', 'C2H2', 'C3H8', 'C3H6', 'C6H6', 'C7H8', 'C10H8']
        plotgasvals=np.array([])
        tar=100-gasdict['N2_MS_avg']-gasdict['H2O_MS_avg']-gasdict['Ar_MS_avg']
        for i in goodgas:
            tar-=gasdict['%s_MS_avg' % i]
            plotgasvals = np.append(plotgasvals, gasdict['%s_MS_avg' % i])
        goodgas.append('C2+')
        plotgasvals = np.append(plotgasvals,tar)

        self.gas_pie_plot = PiePlot(data = plotgasvals, keys = goodgas, figsize = (7,7), save_loc = "%s%s" % (self.directory, str(self.run_info.info['run_id'])))
        self.gas_pie_plot.plot()
        self.gas_pie_plot.save()
        self.gas_pie_plot.close()  
        
        self.run_info.info['piegas']=self.gas_pie_plot.save_loc
        print self.gas_pie_plot.LaTeX_insert('pie_plot_1')
        
    def _load_report_template(self):
        with open('GasificationAnalysisReportTemplate.tex', 'r') as f:
            self.template=f.read()
            
    def generate_standard_report(self):
        self._load_report_template()
        self.text=self.template
        self.variables={i.replace('_', '-'):self.run_info.info[i] for i in self.run_info.info}
        for i in self.variables:
            if type(self.variables[i])==float:
                self.variables[i]='%s' % '%.4g' % self.variables[i]
            else: self.variables[i]=str(self.variables[i])
            self.text=self.text.replace(i, self.variables[i])
        filename=str(self.run_info.info['run_id'])+'_report.tex'
        with open(self.directory+filename,'w') as f:
            f.write(self.text)
        print 'LaTeX file created at %s\n' % self.directory+filename

    def fit_ARIMA(self, col, order = (1,1)):
        try:
            model = ARIMA.ARMA(self.ss[col], order = order)
            result = model.fit()
            self.ss['%s_ARIMA_fitted' % col] = result.fittedvalues
            self.ss['%s_ARIMA_resid' % col] = result.resid
        except KeyError:
            print "Warning: %s is a bad key, ignoring" % col



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = "Run a gasifier analysis")
    parser.add_argument('--run_range', type=str, action = 'store')
    parser.add_argument('--run_id',type=int, action ='store')
    args = parser.parse_args()

    if args.run_id is not None:
        run_id_list = [args.run_id]

    elif args.run_range is not None:
        main_list = args.run_range.split(",")
        run_id_list = []
        for sublist in main_list:
            if ":" in sublist:
                left = sublist.split(":")[0]
                right = sublist.split(":")[1]
                run_id_list.extend(range(int(left), int(right)+1))
            else:
                run_id_list.append(int(sublist))  

    run_id_list=[100]
        #49,50,51,52,53,54,56,57,58,59,
##                 60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,
##                 94,95,96,97,98,99,100,101,102,103,104,105,106,111,112,113,114,120,121,122,123,124,125,126,127,128,129,
##                 130,131,132,133,135,136,137,138,139,140,141,
##                 142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157]
    for run_id in run_id_list:
        print "Generating a Report for run %s..." % run_id
        report = GasifierReport(run_id = run_id)
##    subprocess.call('pdflatex -aux-directory=C:/Users/Ryan.Tracy/Documents/Reports/192 -output-directory=C:/Users/Ryan.Tracy/Desktop C:/Users/Ryan.Tracy/Documents/Reports/192/192_report.tex')


