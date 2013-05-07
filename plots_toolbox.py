import numpy as np
import dataFrame_v2 as df
import matplotlib.pyplot as plt

class Plot:
    """Abstract Data Type for a Plot...this should never be instantiated by itself"""

    def __init__(self, caption = None, figsize = (12,8), save_loc = None, fontsize = 'x-large', subplot = False, subplot_num = 111):
        #Basic plot properties
        self.caption = caption
        self.figsize = figsize
        self.save_loc = save_loc
        self.fontsize = fontsize

        self.subplot = subplot
        self.subplot_num = subplot_num


        #May want to pull in a default dataset here -- we'll see how general this can be
        if subplot:

            plt.subplot(subplot_num)

        else:
            plt.figure(figsize=self.figsize)

    def save(self, loc):
        plt.savefig(loc)

    def close(self):
        plt.close()

    def show(self):
        plt.show()

class PiePlot(Plot):
    """Pie Chart"""

    def __init__(self, data = None, **kwargs):
        Plot.__init__(self,**kwargs)
        if data is None:
            data = OrderedDict()
        if not isinstance(data, OrderedDict):
            raise Exception, "Pie chart values MUST be an ordered dictionary {label:value}"
        self.data = data

    def plot(self):
        
        plt.pie(self.data.items())
        plt.legend(self.data.keys())
        
        plt.tight_layout()

    def save(self):
        loc = self.save_loc+"_gas_comp_pie_plot.png"
        Plot.save(self, loc)

class XYPlot(Plot):
    """This is the basic class for a simple XY plot (one X, many Y)"""
    pass

    def __init__(self, data = None, x_label = None, y_label = None, X_col = None, Y_cols = None, auto_scale = True, legend = True, marker = '-', y_min = None, y_max = None, **kwargs):
        Plot.__init__(self, **kwargs)
        #We should have taken care of whether this is a subplot already...just need to put the plotting machinery in place

        if data is None:
            #data must be a dataframe
            data = df.Dataframe()
        if not isinstance(data, df.Dataframe):
            raise Exception, "XY chart values MUST be in the form of a dataframe"            

        self.data = data

        #More error checking may be appropriate later

        self.X_col = X_col
        self.Y_cols = Y_cols

        if x_label is None and X_col is not None:
            x_label = X_col
        if y_label is None and Y_cols is not None:
            y_label = Y_cols[0]

        self.x_label = x_label
        self.y_label = y_label

        self.auto_scale = auto_scale
        self.legend = legend
        self.marker = marker

        self.y_min = y_min
        self.y_max = y_max

    def plot(self):

        max_val = None
        min_val = None
        legend = []
        plt.xlabel(self.x_label, fontsize = self.fontsize)
        plt.ylabel(self.y_label, fontsize = self.fontsize)
        plt.xticks(fontsize = self.fontsize)
        plt.yticks(fontsize = self.fontsize)

        for y in self.Y_cols:
            legend.append(y)
            
            plt.plot(self.data[self.X_col], self.data[y],self.marker)   #may need a marker default here
            
            
            if max_val is None or np.max(self.data.finite_vals(y)) > max_val:
                max_val = np.max(self.data.finite_vals(y))

            if min_val is None or np.min(self.data.finite_vals(y)) < max_val:
                min_val = np.min(self.data.finite_vals(y))

            if self.auto_scale:
                              
                #set max and min vals to be data members, so that we can get to them later
                self.y_min = 0.97*min_val
                self.y_max = 1.03*max_val

                plt.ylim(self.y_min, self.y_max)
                

            if self.legend:
                plt.legend(legend)
            
            plt.tight_layout()

    def rescale(self, ymin = None, ymax = None, padded=False, padding = 0.03):
        
        if self.subplot:
            plt.subplot(self.subplot_num)
        min_p = 1.0
        max_p = 1.0

        if padded:
            min_p = 1.0-padding
            max_p = 1.0+padding
        if ymin is not None:
            self.y_min = ymin * min_p
        if ymax is not None:
            self.y_max = ymax * min_p
        plt.ylim(self.y_min, self.y_max)

class Histogram(Plot):
    """Creates a histogram for the given data column -- works on a dataframe basis, to allow easy refiguring just by changing column name"""
        
    def __init__(self, data = None, label = None, data_col = None, nbins = 5, useOffset = False, **kwargs):
        Plot.__init__(self,**kwargs)
        if data is None:
            data = df.Dataframe()
        if not isinstance(data, df.Dataframe):
            raise Exception, "data must be a Dataframe!"

        self.data = data
        self.label = label
        self.nbins = nbins
        self.data_col = data_col
        self.useOffset = useOffset

    def plot(self):
        plt.hist(self.data[self.data_col], self.nbins)
        plt.ticklabel_format(useOffset = self.useOffset)
        plt.xlabel(self.label)
        plt.ylabel("Count")

    #inherits save() and close()

class NormalProbabilityPlot(XYPlot):
    """Creates a normal probability plot for the given data colum -- works on a dataframe basis, to allow easy refiguring just by changing column names"""
    def __init__(self, data, data_col = None, **kwargs):
        XYPlot.__init__(self, data = data, **kwargs)
        self.x_label = 'Normal probability plot'
        self.y_label = 'Ordered response'
        self.data_col = data_col
        #Need to set data before this will work

    def _calc_normal_probs(self, data_col):
        n = self.data.numrows()
        U=[1-np.power(0.5,(1/n))]
        ordered=np.sort(self.data[data_col])
        for j in range(0,len(ordered),1)[1:-1]:
            U.append((j-0.3175)/(n+0.365))
        U.append(np.power(0.5,(1/n)))
        U = np.array(U)
        self.data['U_normal_prob'] = U
        self.data['ord_normal_prob'] = ordered

    def plot(self):
        self._calc_normal_probs(self.data_col)
        self.X_col = 'U_normal_prob'
        self.Y_cols = ['ord_normal_prob']
        XYPlot.plot(self)

    def save():
        pass





class LagPlot(XYPlot):
    """Creates a lag plot for the given data column -- works on a dataframe basis, to allow easy refiguring just by changing column names"""
    def __init__(self, data, data_col = None, lag = 1, **kwargs):

        XYPlot.__init__(self, data = data, marker = 'o', **kwargs)
        
        self.x_label = r'Y$_i$'
        self.y_label = r'Y$_{i-%s}$' % lag
        self.data_col = data_col
        self.lag = lag

        
    def _calc_lag(self, data_col):
        for i in range(0, self.lag):
            if i == 0:
                lagged = np.delete(self.data[data_col],0)
            else:
                lagged = np.delete(lagged, 0)
            lagged = np.append(lagged, np.nan)
        self.data['%s_lag' % data_col] = lagged
        
            

    def plot(self):
        self._calc_lag(self.data_col)
        self.X_col = self.data_col
        self.Y_cols = ['%s_lag' % self.data_col]
        XYPlot.plot(self)



class nXYPlot(Plot):
    """Multiple XY subplots on the same plot"""

    def __init__(self, data = None, x_labels = None, y_labels = None, plot_cols = None, h_plots = 1, auto_scale = True, markers = None, **kwargs):
        """Initialize the XY Plot.  data must be a dataframe.  plot_cols is a list of tuples of string/list pairs, with each string corresponding to an X and each list the corresponding Y's to plot"""

        Plot.__init__(self,**kwargs)
        if data is None:
            #Data must be a dataframe
            data = df.Dataframe()
        if not isinstance(data, df.Dataframe):
            raise Exception, "XY chart values MUST be a dataframe"
        self.data = data


        #plot_cols is a list of (string, []), representing (x,Y)
        if plot_cols is None:
            self.plot_cols = []
        self.plot_cols = plot_cols
        
        self.h_plots = h_plots			#For plots with multiple subplots, number of horizontal subplots


        if x_labels is None:
            x_labels = []
            for row in plot_cols:
                x_labels.append(row[0])
            
        if y_labels is None:
            y_labels = []
            for row in plot_cols:
                y_labels.append(row[1][0])

        if not isinstance(x_labels, list):
            raise Exception, "x_labels must be a list of labels, equal to the length of data"
        if not isinstance(y_labels, list):
            raise Exception, "y_labels must be a list of labels, equal to the length of data"
        
        self.markers = markers

       
        self.x_labels = x_labels                
        self.y_labels = y_labels
        

    def plot(self):
        """Plots the X's vs multiple Y's.  Current behavior is to plot each X series in a separate subplot.  Creates separate subplots in a list"""
        plot_list = []
        
        if len(self.plot_cols) % self.h_plots != 0:
            v_plots = len(self.plot_cols)/self.h_plots + 1
        else:
            v_plots = len(self.plot_cols)/self.h_plots
        p_index = 0
        
        #run through and create each of the plots
        for row in self.plot_cols:
            #create a new subplot
          
            subplot_num = v_plots * 100 + self.h_plots * 10 + p_index+1
            
            if self.markers == None:
                marker = '-'
           
            else:
                marker = self.markers[p_index]

            plot_list.append(XYPlot(data = self.data, x_label = self.x_labels[p_index], y_label = self.y_labels[p_index], X_col = row[0], Y_cols = row[1], subplot = True, subplot_num = subplot_num, marker = marker))
            plot_list[p_index].plot()
            p_index += 1
                     
            

class TimeSeriesPlot(nXYPlot):
    """Creates a plot specifically geared to a timeseries"""

    def __init__(self, data = None, Y_cols = None, **kwargs):
        #Y_cols is a list of tuples of the Y vars to be plotted, in the groups indicated by the tuples
        #Put together plot_cols for nXYPlot
        if not isinstance(Y_cols, list):
            raise Exception, "Y_cols must be a list of lists"

        plot_cols = []
        for item in Y_cols:
            if not isinstance(item, list):        
                raise Exception, "Y_cols must be a list of lists"
            for i in item:
                if not isinstance(i, str):
                    raise Exception, "Each item in the lists in Y_cols must be a string"
            plot_cols.append(('timestamp',item))


        nXYPlot.__init__(self, plot_cols = plot_cols, **kwargs)

        if data is not None and not isinstance(data, ts_data):
            raise Exception, "data for a timeseries plot must be an instance of a time series dataframe"

        if not isinstance(plot_cols, list):
            raise Exception, "plot_cols must be a list of columns to plot for a time series plot"
        
        self.data = data

        #self.plot_cols = {'timestamp':plot_cols}

        self.xlabels = ['Time']
        

    #No need to define plot

    def save(self):
        loc = self.save_loc + "_time_series_plot.png"
        Plot.save(self, loc)

    #Need to add fill function and latex interface functions, if necessary

       
        
                    
                
class FourPlot(Plot):
    """Creates a FourPlot from XY data.
       The FourPlot consists of:
       1) A run plot
       2) A lag plot
       3) A histogram
       4) A normal probability plot"""

    #Ideally, I would just build this up from four different plots (two XY, one histogram, ...) but the way I've done subplotting in XY precludes this -- refactoring!

    def __init__(self, data = None, x_label = None, y_label = None, x_var = None, y_var = None, **kwargs):
        #Should put in some x and y variable checking here
        Plot.__init__(self, **kwargs)
        self.x_label = x_label
        self.y_label = y_label
        self.data = data
        self.x_var = x_var
        self.y_var = y_var
        #Need to either store the kwargs or set them for the other plots


    def plot(self):
        """Builds and plots the 4-plot"""

        #Run plot
        self.run_plot = XYPlot(data = self.data, x_label = self.x_label, y_label = self.y_label, X_col = self.x_var, Y_cols = [self.y_var], auto_scale = True, subplot = True, subplot_num = 221)
        self.run_plot.plot()

        #Lag plot
        self.lag_plot = LagPlot(data = self.data, data_col = self.y_var, subplot = True, subplot_num = 222)
        self.lag_plot.plot()

        #Histogram
        self.hist = Histogram(data = self.data, label = self.y_label, data_col = self.y_var, nbins = 20, subplot = True, subplot_num = 223)
        self.hist.plot()

        #Normal probability plot
        self.np_plot = NormalProbabilityPlot(data = self.data, data_col = self.y_var, subplot = True, subplot_num = 224)
        self.np_plot.plot()

    def save(self):
        pass  #Need to implement this



class ControlChartfromDataframe:
    """This is a helper class to take Dataframe data and put it into the appropriate form for a control chart"""

    def __init__(self, data = None, y_col = None, x_col = None, sample_size = 1):
        #Data is in the form of a dataframe
        if data is None:
            data = df.Dataframe()

        if not isinstance(data, df.Dataframe):
            raise Exception, "data must be in the form of a Dataframe"

        self.data = data
        self.y_col = y_col
        self.x_col = x_col
        self.sample_size = sample_size

        

    def getDataframe(self):
            
        grouped_data = [self.data[self.y_col][i-self.sample_size:i] for i in range(0,self.data.numrows(), self.sample_size)[1:]]
        grouped_x = [self.data[self.x_col][i] for i in range(0, self.data.numrows(), self.sample_size)[1:]] # midpoints
        #drop the last group if it is too small -- may want to make this optional
        if len(grouped_data[-1]) != self.sample_size:
            grouped_data.pop()
        self.output_data = df.Dataframe()
        p = 0
        for group, x in zip(grouped_data,grouped_x):
            self.output_data[x] = group
            p += 1
        return self.output_data



class ControlChart(Plot):
    
    #Basic constants
    d2 = [1.128,1.693,2.059,2.326,2.534,2.704,2.847,2.970,3.078,3.173,3.258,3.336,3.407,3.472,3.532,3.588,3.640,3.689,3.735,3.778,3.819,3.858,3.895,3.931]
    d3 = [0.853,0.888,0.880,0.864,0.848,0.833,0.820,0.808,0.797,0.787,0.778,0.770,0.763,0.756,0.750,0.744,0.739,0.734,0.729,0.724,0.720,0.716,0.712,0.708]
    c4 = [0.7979,0.8862,0.9213,0.9400,0.9515,0.9594,0.9650,0.9693,0.9727,0.9754,0.9776,0.9794,0.9810,0.9823,0.9835,0.9845,0.9854,0.9862,0.9869,0.9876,0.9882,0.9887,0.8792,0.9896]



    #All control charts have:
    #    1) Data - a dataframe of the measurements for each sample group
    #    2) Possibly UCL and LCL's
    #    3) Subgroup size? -- really, the common functions are to plot lines given a UCL and LCL, then plot data given the appropriate points
    #    There will always be two graphs, one for x-bar or similar and one for R or s


   
    def __init__(self, data = None, y_label = None, x_label = None, **kwargs):
        Plot.__init__(self, **kwargs)
        self.data = data
        self.y_label = y_label
        self.x_label = x_label
        self.sample_size = self.data.numrows()
        self.chart1_UCL = None
        self.chart1_LCL = None
        self.chart1_target = None

        self.chart2_UCL = None
        self.chart2_LCL = None
        self.chart2_target = None


    def plot(self):
        #set limits and draw horizontal lines
        #plot data on graph from worked up sets
        #common to all kinds of charts

        

        pass        

class XBarControlChart(ControlChart):
        


    def __init__(self, **kwargs):
        ControlChart.__init__(self, **kwargs),

    

    def _calcXBarPoints(self):
        #calculate XBar points
        
        self.X_bar = np.array([])
        self.ord_pts = np.array([])
        for col in self.data.cols():
            self.X_bar = np.append(self.X_bar, self.data[col].mean())   #This is unlikely to work if there are NaN's...need to fix later
            self.ord_pts = np.append(self.ord_pts,col)
        self.X_bar_bar = self.X_bar.mean()

    def plot(self, chart2 = None):
        self._calcXBarPoints()
        self._calcControlLimits()
        data = df.Dataframe(array_dict = {'x-bar':self.X_bar, chart2:getattr(self,chart2), 'x':self.ord_pts})
        self.chart1_plot = XYPlot(data = data, x_label = self.x_label, y_label = 'x-bar', X_col = 'x', Y_cols = ['x-bar'], auto_scale = True, subplot = True, subplot_num = 211, marker = 'o')        
        self.chart1_plot.plot()
        
        if chart2 is not None:
            self.chart2_plot = XYPlot(data = data, x_label = self.x_label, y_label = chart2, X_col = 'x', Y_cols = [chart2], auto_scale = True, subplot = True, subplot_num = 212, marker = 'o')
            self.chart2_plot.plot()


        self.colors = {'UCL':'red', 'LCL':'red', 'target':'blue'}
        for suf in ['UCL', 'LCL', 'target']:
            
            for ch_num in [1,2]:
                plotnum = 210 + ch_num
                plt.subplot(plotnum)
                val = getattr(self, "chart%s_%s" % (ch_num,suf))
                if val is not None:
                    if suf is not 'LCL' or val > 0:
                        plt.hlines(getattr(self, "chart%s_%s" % (ch_num,suf)), min(self.ord_pts), max(self.ord_pts), colors = self.colors[suf])
                        if suf == 'LCL':
                            if getattr(self, "chart%s_LCL" % ch_num) < getattr(self, "chart%s_plot" % ch_num).y_min:
                                getattr(self, "chart%s_plot" % ch_num).rescale(ymin=0.9*getattr(self,"chart%s_LCL" % ch_num), padded = False)
                                
                        elif suf == 'UCL':
                            if getattr(self, "chart%s_UCL" % ch_num) > getattr(self, "chart%s_plot" % ch_num).y_max:
                                getattr(self, "chart%s_plot" % ch_num).rescale(ymax=1.1*getattr(self,"chart%s_UCL" % ch_num), padded = False)

                                    
        

                          
                            
        
  

class XBarRControlChart(XBarControlChart):
    
    def __init__(self, **kwargs):
        XBarControlChart.__init__(self, **kwargs)

    

    def _calcControlLimits(self):
        self.R = np.array([])
        for col in self.ord_pts:
            self.R = np.append(self.R,np.max(self.data[col]) - np.min(self.data[col]))
        self.R_bar = self.R.mean()
        
        #UCL/LCL
        # These could be made adjustable in the future for other than 3-sigma control
        self.chart1_target = self.X_bar_bar
        self.chart1_UCL = self.X_bar_bar + 3./(ControlChart.d2[self.sample_size-2]*np.sqrt(self.sample_size))*self.R_bar
        self.chart1_LCL = self.X_bar_bar - 3./(ControlChart.d2[self.sample_size-2]*np.sqrt(self.sample_size))*self.R_bar

        self.chart2_target = self.R_bar
        self.chart2_UCL = self.R_bar + 3.*ControlChart.d3[self.sample_size-2]/ControlChart.d2[self.sample_size-2]*self.R_bar
        self.chart2_LCL = max(0, self.R_bar - 3.*ControlChart.d3[self.sample_size-2]/ControlChart.d2[self.sample_size-2]*self.R_bar)
        

    def plot(self):
        XBarControlChart.plot(self, chart2 = 'R')
        
        

        

class XBarSControlChart(XBarControlChart):
    
    def __init__(self, **kwargs):
        XBarControlChart.__init__(self, **kwargs)

    def _calcControlLimits(self):

        self.s = np.array([])
        for col in self.ord_pts:
            self.s = np.append(self.s,self.data[col].std())
        self.s_bar = self.s.mean()
        
        #UCL/LCL
        # These could be made adjustable in the future for other than 3-sigma control
        self.chart1_target = self.X_bar_bar
        self.chart1_UCL = self.X_bar_bar + 3./(ControlChart.c4[self.sample_size-2]*np.sqrt(self.sample_size))*self.s_bar
        self.chart1_LCL = self.X_bar_bar - 3./(ControlChart.c4[self.sample_size-2]*np.sqrt(self.sample_size))*self.s_bar

        self.chart2_target = self.s_bar
        self.chart2_UCL = self.s_bar + 3./ControlChart.c4[self.sample_size-2]*np.sqrt(1-np.power(ControlChart.c4[self.sample_size-2],2))*self.s_bar
        self.chart2_LCL = max(0, self.s_bar - 3./ControlChart.c4[self.sample_size-2]*np.sqrt(1-np.power(ControlChart.c4[self.sample_size-2],2))*self.s_bar)

    def plot(self):
        XBarControlChart.plot(self, chart2 = 's')

class IndividualsXBarControlChart(XBarControlChart):

    def __init__(self, **kwargs):
        XBarControlChart.__init__(self, **kwargs)

    def _calcControlLimits(self):
        self.MR = np.array([])
        for index in range(1,len(self.ord_pts)):
            self.MR = np.append(self.MR, np.abs(self.data[self.ord_pts[index]] - self.data[self.ord_pts[index-1]]))
        self.MR_bar = self.MR.mean()
        self.MR = np.insert(self.MR, 0, np.nan) #Need to tack a NaN at the front of the Moving Range, so that it does not plot
        

        #UCL/LCL
        #These could be made adjustable in the future for other than 3-sigma control
        self.chart1_target = self.X_bar_bar
        self.chart1_UCL = self.X_bar_bar + 3.*self.MR_bar/ControlChart.d2[0]
        self.chart1_LCL = self.X_bar_bar - 3.*self.MR_bar/ControlChart.d2[0]

        self.chart2_target = self.MR_bar
        self.chart2_UCL = self.MR_bar + 3.*ControlChart.d3[0]/ControlChart.d2[0]*self.MR_bar
        self.chart2_LCL = max(0, self.MR_bar - 3.*ControlChart.d3[0]/ControlChart.d2[0]*self.MR_bar)

    def plot(self):
        XBarControlChart.plot(self, chart2 = 'MR')
        
