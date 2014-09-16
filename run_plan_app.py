#!/usr/bin/env python

from gi.repository import Gtk
import pandas as pd
import db_toolbox as db
import pandas.io.sql as psql
import os


class PlanApp:

    sp_info = ['temperature', 'pressure', 'biomass_rate', 'steam_flow', 'steam_temp', 'ent_CO2', 'sweep_CO2', 'Ar_tracer', 'superheater_purge', 'tube_diameter']
    bm_info = ['sample_name', 'moisture', 'w_c', 'w_n', 'w_h', 'd10', 'd50', 'd90']
    run_info = ['run_id', 'exp_id', 'ts_start', 'ts_stop', 'ss_start', 'ss_stop', 'operator', 'feeder_slope', 'feeder_intercept']

    labels = {}
    labels['temperature'] = 'Temperature (F)'
    labels['pressure'] = 'Pressure (psig)'
    labels['biomass_rate'] = 'Biomass Flowrate (lb/hr)'
    labels['steam_flow'] = 'Steam Flowrate (lb/hr)'
    labels['steam_temp'] = 'Steam Temperature (F)'
    labels['ent_CO2'] = 'Entrainment CO2 (ACFM)'
    labels['sweep_CO2'] = 'Sweep CO2 (SCFH)'
    labels['Ar_tracer'] = 'Tracer Ar (SCFH)'
    labels['superheater_purge'] = 'Superheater purge (SCFH)'
    labels['tube_diameter'] = 'Tube diameter (in)'
    labels['sample_name'] = 'Biomass name'
    labels['moisture'] = 'Biomass moisture (%)'
    labels['w_c'] = 'wt \% C'
    labels['w_n'] = 'wt \% N'
    labels['w_h'] = 'wt \% H'
    labels['d10'] = 'd10 (micron)'
    labels['d50'] = 'd50 (micron)'
    labels['d90'] = 'd90 (micron)'
    labels['run_id'] = 'run id'
    labels['exp_id'] = 'Exp Tag #'
    labels['ts_start'] = 'Timestamp Start'
    labels['ts_stop'] = 'Timestamp Stop'
    labels['ss_start'] = 'Steady-state Start'
    labels['ss_stop'] = 'Steady-state Stop'
    labels['operator'] = 'Sundrop Supervisor'
    labels['feeder_slope'] = 'Feeder calibration slope'
    labels['feeder_intercept'] = 'Feeder calibration intercept'


    def __init__(self):
        self.builder = Gtk.Builder()
        self.builder.add_from_file("pilot_run_plan_view_GUI.glade")
        self.main_window = self.builder.get_object("window1")
        self.main_window.connect("delete-event", Gtk.main_quit)

        pd = PasswordDialog()
	(rsp, self.user, self.password, self.server) = pd.run()
        
                   

        self._setup_db_connection(self.user, self.password, self.server)
        self._build_treeview()


        self.add_row_button = self.builder.get_object("add_row_button")
        self.add_row_button.connect("clicked", self.add_row)

        self.update_db_button = self.builder.get_object("update_database_button")
        self.update_db_button.connect("clicked", self.update_db)

        self.delete_row_button = self.builder.get_object("delete_row_button")
        self.delete_row_button.connect("clicked", self.delete_row)

        self.run_analysis_button = self.builder.get_object("analyze_data_button")
        self.run_analysis_button.connect("clicked", self.run_analysis)


        self.main_window.show_all()


    def _setup_db_connection(self, user, password, host):
        #set up the database connection
        self.interface = db.db_interface(host = host, user = user, passwd = password)
        self.interface.connect()
        q = db.use_Query("pilot_proc_db")
        self.interface.query(q)


    def _build_treeview(self):
        #get the data into a pandas dataframe
        objects = []
        objects.extend(PlanApp.run_info)
        objects.extend(PlanApp.sp_info)
        objects.extend(PlanApp.bm_info)

        q = db.select_Query(objects = objects, table = "run_plan_view")
        
        self.data = psql.read_frame(q.getQuery(), con = self.interface.connection)


        self.store = Gtk.ListStore(bool, str, str, str, str, str, str, str, str, str, str, str, str, str, str, str, str, str, str, str, str, str, str, str, str, str, str, str) #This is awful -- you really should not do it this way
        #add all of the rows to the store
        self.col_lookup = {}
        i = 1
        for col in objects:
            self.col_lookup[col] = i
            i+=1

        

        for i in range(0, len(self.data.index)):
            #build the row
            row = [False,]
            for col in objects:
                row.append(str(self.data.iloc[i][col]))
            
            self.iter = self.store.append(row)

            self.row_count = i


        #now we want to build the view
        self.view = self.builder.get_object("run_plan_table")
        self.view.set_model(self.store)

        self.view_renderers = {}
        self.view_columns = {}
        
        self.analysis_renderer = Gtk.CellRendererToggle()
        self.analysis_col = Gtk.TreeViewColumn("Analyze", self.analysis_renderer, active = 0)
        self.analysis_renderer.connect("toggled", self.on_analysis_toggled)
        self.analysis_renderer.set_sensitive(True)
        self.view.append_column(self.analysis_col)

        for col in objects:
            self.view_renderers[col] = Gtk.CellRendererText()
            self.view_renderers[col].set_property("editable", True)
            self.view_renderers[col].connect("edited", self.text_edited, self.col_lookup[col])


            self.view_columns[col] = Gtk.TreeViewColumn(PlanApp.labels[col], self.view_renderers[col], text = self.col_lookup[col])
            self.view.append_column(self.view_columns[col])
         
        self.select = self.view.get_selection()
        self.select.connect("changed", self.on_tree_selection_changed)
           
    def on_analysis_toggled(self, widget, path):
        self.store[path][0] = not self.store[path][0]

    def text_edited(self, widget, path, text, col_num = 0):
        self.store[path][col_num] = text
        
    def on_tree_selection_changed(self, selection):
        
        model, self.current_iter = self.select.get_selected()
        

    def add_row(self, widget):
        self.iter = self.store.append()
        self.row_count += 1
        self.store[self.iter][1] = str(self.row_count) 

    def delete_row(self, path):
        #print self.store[self.current_iter]
        if self.current_iter is not None:
            self.store.remove(self.current_iter)
    
    def run_analysis(self, widget):
        #iterate through the model and get the run ids
        run_ids = []
        it = self.store.get_iter_first()
        ok = self.store.iter_is_valid(it)
        if ok:

            while ok:
                if self.store[it][0] == True:
                    run_ids.append(self.store[it][1])
                it = self.store.iter_next(it)
                ok = (it is not None)

        run_list = ",".join(run_ids)
        os.system("python ./pilot_analysis.py --run_range %s --user %s --host %s --pswd %s" % (run_list, self.user, self.server, self.password))
                       


    def update_db(self, widget):
        self._generate_CSV_file()
        os.system('python ./run_plan_fill.py --filename proxy_csv.csv --host %s --user %s --pswd %s' % (self.server, self.user, self.password))
        print "done updating"

    def _generate_CSV_file(self):
        f = open('proxy_csv.csv', 'w')
        #get the initial iterator
        header = ",".join(PlanApp.run_info)
        header += ',' + ",".join(PlanApp.sp_info)
        header += ',' + ','.join(PlanApp.bm_info)
        header += "\n"
        f.write(header)
        
        
        it = self.store.get_iter_first()
        ok = self.store.iter_is_valid(it)
        if ok:
            
            while ok:
                row = ""
                rl = []
                for item in self.store[it]:
                    rl.append(item)
                rl = rl[1:]
                for item in rl:

                    if item is not None:
                        row += "%s," % item
                    else:
                        row += "NULL,"
                row = row[:-1] + "\n"
                f.write(row)
                it = self.store.iter_next(it)
                ok = (it is not None)
                       

        f.close()


        print "Done writing csv"

    def _call_fill_code(self):
        pass


    def run(self):
        Gtk.main()


class PasswordDialog:

    def __init__(self):
        self.builder = Gtk.Builder()
        self.builder.add_from_file('mysql_login_dialog.glade')
        self.dialog_window = self.builder.get_object('dialog1')

        self.login_button = self.builder.get_object('login_button')
        self.cancel_button = self.builder.get_object('cancel_button')

        self.login_button.connect("clicked", self.login)
        self.cancel_button.connect("clicked", self.cancel)

        #self.dialog_window.run()


    def run(self):
        self.dialog_window.run()
        return (self.response, self.login, self.password, self.host)

    def login(self, widget):
        le = self.builder.get_object('login_entry')
        pe = self.builder.get_object('pass_entry')
        he = self.builder.get_object('host_entry')
        self.login = le.get_text()
        self.password = pe.get_text()
        self.host = he.get_text()
        if self.host == '' or self.host is None:
            self.host = 'localhost'
        self.response = True

        self.dialog_window.destroy()

    def cancel(self, widget):
        print "cancel clicked"
        self.login = None
        self.password = None
        self.host = None
        self.response = False
        self.dialog_window.destroy()


if __name__ == "__main__":

    app = PlanApp()
    app.run()
