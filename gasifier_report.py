#!/usr/bin/env python

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
        self.interface_raw = db.db_interface(host = "192.168.13.51", user = "dbmaint", passwd = "f9p2#nH1")
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
   
        
        #Generate pie plot
        
        self.calc_gas_comp_pie_plot()
        self.gas_pie_plot = PiePlot(data = self.pie_gasvals, keys = self.pie_goodgas, figsize = (7,7), save_loc = "%s%s" % (self.directory, str(self.run_info.info['run_id'])))
        self.gas_pie_plot.plot()
        self.gas_pie_plot.save()
        pie_LaTeX = self.gas_pie_plot.save_loc
        self.gas_pie_plot.close()
          
        
        self.run_info.info['piegas']=pie_LaTeX


        ts_plots = {}
        #Generate time series plots
        ts_keys = ['mass_feed','main_comp','trace_comp','tube_temps']
        ts_Ycols = [['mass_flow_brush_feeder'],['CO_MS','CO2_MS', 'H2_MS', 'CH4_MS','Ar_MS'],['C2H4_MS', 'C6H6_MS', 'C7H8_MS', 'C10H8_MS'],['temp_skin_tube_top','temp_skin_tube_middle','temp_skin_tube_bottom']]
        ts_ylabels = ['Biomass feed rate (lb/hr)', 'Gas Composition (% vol)','Gas Composition (ppm)', 'Tube Skin Temperatures ($^\circ$C)']
        ts_captions = ['Time series plot for biomass flowrate', 'Time series plot for gas composition', 'Time series plot for gas composition', 'Time series plot for reactor tube skin temperatures']
        ts_markers = ['-','o','o','-']
        
        LaTeX_ts = ""

        for key, cols, label, caption, marker in zip(ts_keys,ts_Ycols, ts_ylabels, ts_captions, ts_markers):
            ts_plots[key] = TimeSeriesPlot(data = self.ts, Y_cols = [cols], y_labels = [label], caption = caption, save_loc = "%s%s_%s" % (self.directory, self.run_id, key), markers =[marker])
            ts_plots[key].plot()
            ts_plots[key].fill(self.ss['timestamp'])
            ts_plots[key].save()
            LaTeX_ts += ts_plots[key].LaTeX_insert("ts_%s" % key)
            ts_plots[key].close()

        #Generate four plots
        fp_plots = {}
        fp_keys = ['mass_feed', 'temp_mid', 'temp_steam','pressure_KO', 'CO_MS', 'CO2_MS', 'H2_MS', 'CH4_MS']
        fp_Y = ['mass_flow_brush_feeder','temp_skin_tube_middle','temp_steam_reactor_entry','pressure_ash_knockout_vessel','CO_MS','H2_MS','CO2_MS','CH4_MS']
        fp_label = ['Biomass Flow Rate (lbs/hr)','Reactor Skin Middle ($^\circ$C)','Steam Temperature ($^\circ$C)','Ash Knockout Pressure (psig)','Carbon Monoxide (mol%)','Hydrogen (mol%)','Carbon Dioxide (mol%)', 'Methane (mol%)']
        fp_caption = ['Four-plot for biomass flow rate','Four-plot for reactor skin temperature','Four-plot for temperature of steam at reactor inlet', 'Four-plot for ash knockout pressure', 'Four-plot for carbon monoxide readings from the mass spectrometer','Four-plot for hydrogen readings from the mass spectrometer','Four-plot for carbon dioxide readings from the mass spectrometer','Four-plot for methane readings from the mass spectrometer']
        
        LaTeX_fp = ""

        for key, Y, label, caption in zip(fp_keys, fp_Y, fp_label, fp_caption):
            fp_plots[key] = FourPlot(data = self.ss, x_label = 'Time', y_label = label, x_var = 'timestamp', y_var = Y, caption = caption,save_loc = "%s%s_%s" % (self.directory, self.run_id, key))
            fp_plots[key].plot()
            fp_plots[key].save()
            LaTeX_fp += fp_plots[key].LaTeX_insert("fp_%s" % key)
            fp_plots[key].close()
        
                
        #Create ARIMA fits as necessary - will not work for NaN data (i.e. raw MS data -- that should not be autocorrellated anyway)
        ARIMA_list = ['mass_flow_brush_feeder','temp_steam_reactor_entry']
        for col in ARIMA_list:
            #try:
                self.fit_ARIMA(col)
                
            #except Exception:
            #    ARIMA_list.remove(col)

        ARIMA_captions = {'mass_flow_brush_feeder':'biomass flow rate','temp_steam_reactor_entry':'reactor inlet steam temperature'}        

        cp_plots = {}
        cp_keys = ['mass_feed', 'temp_mid', 'temp_steam','pressure_KO', 'CO_MS', 'CO2_MS', 'H2_MS', 'CH4_MS', 'mass_feed_ARIMA','temp_steam_ARIMA']
        cp_Y = ['mass_flow_brush_feeder','temp_skin_tube_middle','temp_steam_reactor_entry','pressure_ash_knockout_vessel','CO_MS','H2_MS','CO2_MS','CH4_MS']
        cp_caption = []
        items = ['biomass flow rate', 'reactor skin temperature', 'temperature of steam at reactor inlet', 'ash knockout pressure', 'carbon monoxide (MS)', 'hydrogen (MS)' ,'carbon dioxide (MS)', 'methane (MS)']
        for col in ARIMA_list:
            cp_keys.append('%s_ARIMA' % col)
            cp_Y.append('%s_ARIMA_resid' % col)
            items.append(ARIMA_captions[col])
        for item in items:
            cp_caption.append("Individuals control chart for %s ARIMA(1,1) residuals" % item)
        
            


        LaTeX_cp = ""
        for key, Y, caption in zip(cp_keys, cp_Y, cp_caption):
            input_df = ControlChartfromDataframe(data = self.ss, y_col = Y, x_col = 'timestamp', ignore_nan = True)
            cp_plots[key] = IndividualsXBarControlChart(data = input_df.getDataframe(), caption = caption, save_loc = "%s%s_%s" % (self.directory, self.run_id, key))
            cp_plots[key].plot()
            cp_plots[key].annotate(1)
            cp_plots[key].annotate(2)
            cp_plots[key].save()
            LaTeX_cp += cp_plots[key].LaTeX_insert("cp_%s" % key)
            cp_plots[key].close()

        self.run_info.info['timeseriesplots'] =  LaTeX_ts
        self.run_info.info['fourplots'] = LaTeX_fp
        self.run_info.info['controlplots'] = LaTeX_cp
           
        
       

        
        self.generate_standard_report()
        
        
        
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
        
    
      
    def calc_gas_comp_pie_plot(self):
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
        self.pie_gasvals = np.append(plotgasvals,tar)
        self.pie_goodgas = goodgas
        
        
        
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
        with open('./'+filename,'w') as f:
            f.write(self.text)
        print 'LaTeX file created at %s\n' % '/.'+filename

    def fit_ARIMA(self, col, order = (1,1)):
        try:
            
            model = ARIMA.ARMA(self.ss[col])
            
            result = model.fit(order=order)
            self.ss['%s_ARIMA_fitted' % col] = result.fittedvalues
            self.ss['%s_ARIMA_resid' % col] = result.resid
        except KeyError:
            print "Warning: %s is a bad key, ignoring" % col

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
        print "Generating a Report for Run %s..." % run_id
        try: report = GasifierReport(run_id = run_id)
        except: print "Report Generation Failed for Run %s" %run_id



