import LabFunctionLib as lfl
import numpy as np
import db_toolbox as SQL
import db_toolbox as db
import datetime

import pygtk
pygtk.require('2.0')
import gtk

class Lab_Gasifier_GUI:
    def open_grapher(self, widget, data=None):
        df=self.generate_dataFrame()
        Graph_GUI(df)
    
    def generate_dataFrame(self):
        start = self.startentry.get_text()
        end=self.endentry.get_text()
        
        dbconn = db.db_interface(host = "192.168.13.51", user = "dbmaint", passwd = "f9p2#nH1")
        dbconn.connect()
        q = db.use_Query("lab_run_db")
        dbconn.query(q)
        
        dtstart=datetime.datetime.strptime(start, '%Y-%m-%d %H:%M:%S')
        dtend=datetime.datetime.strptime(end, '%Y-%m-%d %H:%M:%S')
        df=lfl.ts_data(dtstart, dtend)
        df.SQL_load(dbconn, table='gasifier_data_HR_view')
##        df.get_ms_from_csv('12101214.csv')
        df.interp_ga()
##        df.calc_outlet_flowrate()
##        df.calc_gas_produced_flowrate()
##        df.calc_product_flowrates()
##        df.calc_carbon_in(0.522)
##        df.calc_carbon_out()
        return df
        
    def __init__(self):
        window = gtk.Window(gtk.WINDOW_TOPLEVEL)
        window.set_title("Laboratory Gasification")
        window.connect("delete_event", lambda w,e: gtk.main_quit())
        
        table=gtk.Table(rows=3, columns=2)
        window.add(table)

        startlabel=gtk.Label('Start')
        table.attach(startlabel, 0, 1, 0, 1)
        startlabel.show()

        endlabel=gtk.Label('End')
        table.attach(endlabel, 0, 1, 1, 2)
        endlabel.show()

        self.startentry=gtk.Entry()
        self.startentry.set_text('2012-12-10 15:50:00')
        table.attach(self.startentry, 1, 2, 0, 1)
        self.startentry.show()

        self.endentry=gtk.Entry()
        self.endentry.set_text('2012-12-10 16:00:00')
        table.attach(self.endentry, 1, 2, 1, 2)
        self.endentry.show()

        graphbutton=gtk.Button('Graph')
        graphbutton.connect('clicked', self.open_grapher)
        table.attach(graphbutton, 0, 2, 2, 3)
        graphbutton.show()
        
        table.show()
        window.show()
        
    
class Graph_GUI:
    def check(self, widget, data=None):
        if widget.get_active()==True:
            self.graphlist.append(data)
        if widget.get_active()==False:
            del(self.graphlist[self.graphlist.index(data)])

    def graph(self, widget, data=None):
        import matplotlib.pyplot as plt
        for i in self.graphlist:
            plt.plot(data['timestamp'], data[i], label=i.replace('_', ' '))
        plt.figure(1).autofmt_xdate()
        plt.title(self.graphtitle.get_text())
        plt.ylabel(self.graphylabel.get_text())
        plt.xlabel('Time')
        plt.legend(loc='best')
        plt.show()
        
    def __init__(self, df):
        keydict=df.list_keys()
        keylist=keydict.keys()
        self.graphlist=[]
        ncols=6
        cols=range(ncols)
        nrows=(len(keydict)/ncols+1)*2+3
        oddrows=[]
        evenrows=[]
        buttonlist={}

        graphwindow=gtk.Window(gtk.WINDOW_TOPLEVEL)
        graphwindow.set_title('Laboratory Gasication')
        graphwindow.set_border_width(10)
        
        table=gtk.Table(rows=nrows, columns=ncols)
        graphwindow.add(table)
        
        for i in range(nrows):
            if round(i/2)==i/2.:
                evenrows.append(i)
            else:
                oddrows.append(i)
        for i in range(len(keylist)+len(evenrows)):
            try:
                label=gtk.Label('\n'+keylist[0]+'\n')
                try:
                    colstart=cols[0]
                    colend=cols[0]+1
                except IndexError:
                    pass
                rowstart=evenrows[0]
                rowend=evenrows[0]+1
                try:
                    del(cols[0])
                    table.attach(label, colstart, colend,
                                 rowstart, rowend)
                    label.show()
                    vbox=gtk.VBox(True, 0)
                    for b in keydict[keylist[0]]:
                        buttonlist[b]=gtk.CheckButton(b.replace('_',' '))
                        buttonlist[b].connect("toggled", self.check, b)
                        vbox.pack_start(buttonlist[b], True, True, 0)
                        buttonlist[b].show()
                    table.attach(vbox, colstart, colend, rowstart+1, rowend+1)
                    vbox.show()
                    del(keylist[0])                                                  
                except IndexError:
                    del(evenrows[0])
                    cols=range(ncols)
            except IndexError: pass

        label=gtk.Label('\n Graph Title \n')
        table.attach(label, 0, 1, nrows-3, nrows-2)
        label.show()
        self.graphtitle=gtk.Entry()
        table.attach(self.graphtitle, 1, 2, nrows-3, nrows-2)
        self.graphtitle.show()

        label=gtk.Label('\n Y-axis Label \n')
        table.attach(label, 0, 1, nrows-2, nrows-1)
        label.show()
        self.graphylabel=gtk.Entry()
        table.attach(self.graphylabel, 1, 2, nrows-2, nrows-1)
        self.graphylabel.show()
        
        graphbutton=gtk.Button('Create Graph')
        graphbutton.connect('clicked', self.graph, df)
        table.attach(graphbutton, 0, 2, nrows-1, nrows)
        graphbutton.show()
        table.show()
        graphwindow.show()
        
    
def main():
    gtk.main()
    return 0       

if __name__ == "__main__":
    Lab_Gasifier_GUI()
    main()
