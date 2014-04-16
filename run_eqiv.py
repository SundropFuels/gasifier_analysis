import numpy as np
import dataFrame_pd as df
import db_toolbox as db
import getpass
import scipy.interpolate as spi
import scipy.stats as st

class eq_set:

    #Use of the dictionaries may be a bit slow -- could use dataframes if we did not need to dynamically allocate or if we left a lot of space
    """
    def __init__(self, base_set = None):
        if base_set is None:
            base_set = {}
        
        if not isinstance(base_set, dict) and not isinstance(base_set, tuple):
            raise TypeError, "The base set must be a dictionary or a tuple"
        if isinstance(base_set, tuple):
            base_set = {base_set[0]: base_set[1]}
        self.points = base_set  #This contains labels with values for categorization
        
        if len(base_set) != 0:
            self.label = np.array(base_set.values()).mean()
    """
    def __init__(self, base_set = None):
        
        if not isinstance(base_set, dict):
            raise TypeError, "The base set MUST be a dictionary with an 'id' field and a 'value' field"
        #implementation of this now is with a pair of numpy arrays -- not necessarily fast, but the way the finder works we will add whole sets of lists               
        self.points = base_set
        if len(self.points['value']) != 0:
            self.label = self.points['value'].mean()

    """
    def extend(self, point):
        if not isinstance(point, tuple):
            raise TypeError, "The data point needs to be a tuple"
        self.points[point[0]]=point[1]
        #need to update the label to include the new point
        self.label = np.array(self.points.values()).mean()
        print "added run# %s.\tIt's val: %s.\tNew label: %s" % (point[0], point[1], self.label)
    """
    def extend(self, points):
        if not isinstance(points, tuple):
            raise TypeError, "The data points need to be a tuple of numpy arrays"
        for item in points:
            if not isinstance(item, np.array):
                raise TypeError, "The data points and ids need to be numpy arrays"
        self.points['id'] = np.append(self.points['id'], points[0])
        self.points['value'] = np.append(self.points['value'], points[1])
        self.label = self.points['value'].mean()
    """
    def contains(self, point, method = "reldiff", rtol=0.1, atol=0.1):
        if not isinstance(point, tuple):
            raise TypeError, "The data point needs to be a tuple"
        if self.label != 0:
            tol = rtol
        else:
            tol = atol
        return getattr(self, "_%s" % method)(point[1]) < tol

    def _reldiff(self, val):
        #print self.label
        if self.label != 0:
            return abs((val - self.label))/self.label
        else:
            return abs(val - self.label)

    def _reldiff_corr(self, val, corr):
        if self.label != 0;
            return abs((val - self.label))
    """
    def output(self):
        print "Equivalence set: %s" % self.label
        for key in self.points:
            print "%s:\t%s" % (key, self.points[key])

    def keys(self):
        return self.points['id'].tolist()

class partitionDataframe(df.Dataframe):

    def __init__(self, **kwargs):
        df.Dataframe.__init__(self, **kwargs)
        self.partition_tolerance = 0.1
        self.partition_tolerance_abs = 0.1

    def partition(self, cols, id_col):
        if not isinstance(cols, list):
            raise TypeError, "The list of columns to partition on must be a list object" ##!!FINISH 

        eq_sets = {}
        
        for col in cols:
            """
            i = 1
            #Need to walk the data and determine which equivalence set each piece belongs to
            eq_sets[col] = [eq_set(base_set = (self[id_col][0], self[col][0]))]
            #We will want to subtract away the mean value for all of the data sets so that the relative tolerance better represents deviation
            mean = self[col].mean()
            while i < len(self.index):
                #check if the next point is in the existing list of equivalence sets -- stop when you find one
                found = False
                for es in eq_sets[col]:
                    point = (self[id_col][i], self[col][i])
                    if es.contains(point, rtol = self.partition_tolerance, atol = self.partition_tolerance_abs):
                        found = True
                        es.extend(point)
                        break	#quit early -- we found the right enclosing set
                if not found:
                    eq_sets[col].append(eq_set(point))  #add a new equivalence set to the list
                i += 1
            """
            eq_sets[col] = self._create_equivalence_sets(col, id_col)
        self.eq_sets = eq_sets

    def _create_equivalence_sets(self, col, id_col):
        """Creates an equivalence set using a kernel density method for the given column"""
        #Create the kernel density estimate.  For starters, I am going to use the default bandwidth (although it may not work as well as desired)
        
        k = st.gaussian_kde(self[col], bw_method = 0.1)
        #Create a set from which we can find the maxima and minima
        h = (self[col].max() - self[col].min())/10000.0 #want 10000 data points in the interpolated set
        x = np.arange(self[col].min(), self[col].max(), h)
	#Estimate the derivatives
        x_plus = x+h
        x_minus = x-h
        
        
        dk = (k(x_plus) - k(x_minus))/(2.0*h)
        d2k = (k(x_plus) - 2.0*k(x) + k(x_minus))/(h**2)
        #Fit splines, find minima
        dks = spi.InterpolatedUnivariateSpline(x, dk)
        r = dks.roots()
        d2ks = spi.InterpolatedUnivariateSpline(x,d2k)
        bounds = [self[col].min(),]
        
        for root in r:
            if d2ks(root)>0:		#second derivative > 0 means a minimum
                bounds.append(root)
            
        bounds.append(self[col].max()+1.0)
        intervals = {}
        print bounds
        i = 0
        eq_sets = []
        while i < len(bounds)-1:
            
            #find the points in each interval
            
             
            ids = self[id_col][np.logical_and(self[col]>= bounds[i],self[col]<bounds[i+1])]
            vals = self[col][np.logical_and(self[col]>= bounds[i],self[col]<bounds[i+1])]
            
            es_new = eq_set(base_set = {'id':ids, 'value':vals})
            
            eq_sets.append(es_new)
            i += 1
        return eq_sets

        
        


    def unique_set(self, lists):
        """Finds the set of points with unique conditions on the partition"""
        #We want this to return a list of lists of run numbers
        
	#base case
        if len(lists) == 1:
            
            return lists[0]

        #partition out the first item
        
        listsA = lists[0]
        listsB = self.unique_set(lists[1:])

        #walk through the lists to create the unique set of items
        new_lists = []
        for listA in listsA:
            for listB in listsB:
                
                new_list = []
                for itemA in listA:
                    
                    if itemA in listB:
                        new_list.append(itemA)
                        #listA.remove(itemA)          #should drop itemA from listA so that we don't look for it twice
                if new_list:			#if the new list is not empty, add it to the list of new lists
                    new_lists.append(new_list)
        return new_lists

    def find_unique_sets(self, cols, id_col):
        self.partition(cols, id_col)
        #Create the list of lists
        l = []
        for col in self.eq_sets:
            li = []
            for es in self.eq_sets[col]:
                li.append(es.keys())
            l.append(li)
        return self.unique_set(l)


class EquivalentSetFinder:

    def __init__(self, user, password):
        #Create the interface to the data and load the partition data frame
        self.interface = db.db_interface(host = "192.168.13.51", user = user, passwd = password)
        self.interface.connect()
        q = db.use_Query("lab_proc_db")
        self.interface.query(q)

        self.runs = partitionDataframe()
        self.runs.SQL_load_data(self.interface, table = 'integral_summary')

        
    def find_unique_sets(self, cols):
        #create the partition of the data frame
        self.unique = self.runs.find_unique_sets(cols, "run_id")
        #for es in self.runs.eq_sets.values():
        #    print es
        #    for i in range(0,len(es)):
        #        es[i].output()

        #due to ordering considerations, these lists may not be 

        #build the list of lists
        for un in self.unique:
            print un
           
            

        
if __name__ == "__main__":
    #pdf = partitionDataframe()
    #lists = [[[2,3,4],[6,7],[5,8,4],[1]],[[2,3],[6,4],[5,7,8],[1]],[[1,2,3],[4,5,6],[7,8]]]
    #print pdf.unique_set(lists)
    user = raw_input('User: ')
    pswd = getpass.getpass()

    finder = EquivalentSetFinder(user, pswd)
    finder.find_unique_sets(cols = ["space_time_avg"])#, "d50", "pp_H2O_avg", "tube_dia", "pp_CO2_avg", "pressure_ash_knockout_vessel_avg", "mass_flow_brush_feeder_avg","temp_skin_tube_middle_avg"])