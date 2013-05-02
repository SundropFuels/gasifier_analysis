import LabFunctionLib as lfl
import numpy as np
import matplotlib.pyplot as plt
import db_toolbox as SQL
import datetime

import pygtk
pygtk.require('2.0')
import gtk


class Table:
    def delete_event(self, widget, event, data=None):
        gtk.main_quit()
        return False

    def close_event(self, widget, event, data=None):
        self.graphwindow.destroy()
        return False
    
    def check(self, widget, data=None):
        if widget.get_active()==True:
            self.graphlist.append(data)
        if widget.get_active()==False:
            del(self.graphlist[self.graphlist.index(data)])
        print self.graphlist

    def graph(self, widget, data=None):
        for i in self.graphlist:
            plt.plot(self.data['timestamp'], self.data[i])
        plt.show()
        
        
    def callback(self, widget):
        start = self.startentry.get_text()
        end=self.endentry.get_text()
        print start
        print end
        dbconn=SQL.db_interface(host='io', user='jmp_user', passwd='jmpyourface')
        dbconn.connect()
        dtstart=datetime.datetime.strptime(start, '%Y-%m-%d %H:%M:%S')
        dtend=datetime.datetime.strptime(end, '%Y-%m-%d %H:%M:%S')
        self.data=lfl.ts_data(dtstart, dtend)
        self.data.SQL_load(dbconn, table='gasifier_data_HR_view')
        self.keylist=self.data.list_keys()
        self.graphwindow = gtk.Window(gtk.WINDOW_TOPLEVEL)
        self.graphwindow.set_title("Laboratory Gasification")
        self.graphwindow.connect("delete_event", self.close_event)
        hbox=gtk.HBox(True, 0)
        self.graphwindow.add(hbox)
        self.graphwindow.show()
        hbox.show()
        self.buttonlist={}
        self.graphlist=[]
        vboxlist={}
        for i in self.keylist:
            vboxlist[i]=gtk.VBox(False, 0)
            hbox.add(vboxlist[i])
            
            label=gtk.Label(i)
            vboxlist[i].pack_start(label, True, True, 0)
            label.show()
            for b in self.keylist[i]:
                self.buttonlist[b]=gtk.CheckButton(b.replace('_',' '))
                self.buttonlist[b].connect("toggled", self.check, b)
                vboxlist[i].pack_start(self.buttonlist[b], True, True, 0)
                vboxlist[i].show()
                self.buttonlist[b].show()
        graphbutton=gtk.Button('Create Graph')
        graphbutton.connect("clicked", self.graph, self.graphlist)
        self.graphwindow.add(hbox)
        hbox.show()
        hbox.pack_start(graphbutton, False, False, 0)
        graphbutton.show()
        
            

    def __init__(self):
        # create a new window
        window = gtk.Window(gtk.WINDOW_TOPLEVEL)
        window.set_title("Laboratory Gasification")
        window.connect("delete_event", lambda w,e: gtk.main_quit())

        self.vbox = gtk.VBox(True, 0)
        window.add(self.vbox)
        self.vbox.show()

        hbox1=gtk.HBox(False, 0)
        self.vbox.add(hbox1)
        hbox1.show()
        hbox2=gtk.HBox(False, 0)
        self.vbox.add(hbox2)
        hbox2.show()
        
        startlabel = gtk.Label('Start')
        startlabel.set_justify(gtk.JUSTIFY_RIGHT)
        hbox1.pack_start(startlabel, True, True, 0)
        startlabel.show()
        self.startentry = gtk.Entry()
        self.startentry.set_text("2012-12-10 16:05:00")
        hbox1.pack_start(self.startentry, True, True, 0)
        self.startentry.show()

        endlabel = gtk.Label('End ')
        hbox2.pack_start(endlabel, True, True, 0)
        endlabel.show()
        endlabel.set_justify(gtk.JUSTIFY_RIGHT)
        self.endentry = gtk.Entry()
        self.endentry.set_text("2012-12-10 16:10:00")
        hbox2.pack_start(self.endentry, True, True, 0)
        self.endentry.show()

        hbox3 = gtk.HBox(False, 0)
        self.vbox.add(hbox3)
        hbox3.show()
                                                            
        button = gtk.Button("Graph")
        button.connect("clicked", self.callback)
        self.vbox.pack_start(button, True, True, 0)
        button.set_flags(gtk.CAN_DEFAULT)
        button.grab_default()
        button.show()
        window.show()


class CheckButton:
    # Our callback.
    # The data passed to this method is printed to stdout
    def callback(self, widget, data=None):
        print "%s was toggled %s" % (data, ("OFF", "ON")[widget.get_active()])

    # This callback quits the program
    def delete_event(self, widget, event, data=None):
        gtk.main_quit()
        return False

    def __init__(self):
        # Create a new window
        self.window = gtk.Window(gtk.WINDOW_TOPLEVEL)

        # Set the window title
        self.window.set_title("Check Button")

        # Set a handler for delete_event that immediately
        # exits GTK.
        self.window.connect("delete_event", self.delete_event)

        # Sets the border width of the window.
        self.window.set_border_width(20)

        # Create a vertical box
        vbox = gtk.VBox(True, 2)

        # Put the vbox in the main window
        self.window.add(vbox)

        # Create first button
        button = gtk.CheckButton("check button 1")

        buttonlist={'check1':gtk.CheckButton('check1'),
                    'check2':gtk.CheckButton('check2'),
                    'check3':gtk.CheckButton('check3'),
                    'check4':gtk.CheckButton('check4')}

        # When the button is toggled, we call the "callback" method
        # with a pointer to "button" as its argument
        button.connect("toggled", self.callback, "check button 1")


        # Insert button 1
        vbox.pack_start(button, True, True, 2)

        button.show()

        # Create second button

        button = gtk.CheckButton("check button 2")

        # When the button is toggled, we call the "callback" method
        # with a pointer to "button 2" as its argument
        button.connect("toggled", self.callback, "check button 2")
        # Insert button 2
        vbox.pack_start(button, True, True, 2)

        button.show()

        for b in buttonlist:
            buttonlist[b].connect("toggled", self.callback, b)
            vbox.pack_start(buttonlist[b], True, True, 2)
            buttonlist[b].show()

        # Create "Quit" button
        button = gtk.Button("Quit")

        # When the button is clicked, we call the mainquit function
        # and the program exits
        button.connect("clicked", lambda wid: gtk.main_quit())

        # Insert the quit button
        vbox.pack_start(button, True, True, 2)

        button.show()
        vbox.show()
        self.window.show()

def main():
    gtk.main()
    return 0       

if __name__ == "__main__":
    Table()
    main()
