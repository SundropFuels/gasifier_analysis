import numpy as np
import dataFrame_pd as df
import db_toolbox as db
import getpass


class eq_set:

    #Use of the dictionaries may be a bit slow -- could use dataframes if we did not need to dynamically allocate or if we left a lot of space

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
             
        
    def extend(self, point):
        if not isinstance(point, tuple):
            raise TypeError, "The data point needs to be a tuple"
        self.points[point[0]]=point[1]


    def contains(self, point, method = "reldiff", tol=0.1):
        if not isinstance(point, tuple):
            raise TypeError, "The data point needs to be a tuple"
        return getattr(self, "_%s" % method)(point[1]) < tol

    def _reldiff(self, val):
        if self.label != 0:
            return (val - self.label)/self.label
        else: return 0

    def output(self):
        print "Equivalence set: %s" % self.label
        for key in self.points:
            print "%s:\t%s" % (key, self.points[key])

class partitionDataframe(df.Dataframe):

    def __init__(self, **kwargs):
        df.Dataframe.__init__(self, **kwargs)
        self.partition_tolerance = 0.1

    def partition(self, cols):
        if not isinstance(cols, list):
            raise TypeError, "The list of columns to partition on must be a list object" ##!!FINISH 

        eq_sets = {}
        
        for col in cols:
            i = 1
            #Need to walk the data and determine which equivalence set each piece belongs to
            eq_sets[col] = [eq_set(base_set = (self.index[0], self[col][0]))]
            while i < len(self.index):
                #check if the next point is in the existing list of equivalence sets -- stop when you find one
                found = False
                for es in eq_sets[col]:
                    point = (self.index[i], self[col][i])
                    if es.contains(point):
                        found = True
                        es.extend(point)
                        break	#quit early -- we found the right enclosing set
                if not found:
                    eq_sets[col].append(eq_set(point))  #add a new equivalence set to the list
                i += 1

        self.eq_sets = eq_sets

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

class EquivalentSetFinder:

    def __init__(self, user, password):
        #Create the interface to the data and load the partition data frame
        self.interface = db.db_interface(host = "192.168.13.51", user = user, passwd = password)
        self.interface.connect()
        q = db.use_Query("lab_proc_db")
        self.interface.query(q)

        self.runs = partitionDataframe()
        self.runs.SQL_load_data(self.interface, table = 'gas_integral_tbl')

        #print self.runs.index
        #print self.runs.columns

    def find_unique_sets(self, cols):
        #create the partition of the data frame
        self.runs.partition(cols)
        for es in self.runs.eq_sets.values():
            es[0].output()
            
            

        
if __name__ == "__main__":
    #pdf = partitionDataframe()
    #lists = [[[2,3,4],[6,7],[5,8,4],[1]],[[2,3],[6,4],[5,7,8],[1]],[[1,2,3],[4,5,6],[7,8]]]
    #print pdf.unique_set(lists)
    user = raw_input('User: ')
    pswd = getpass.getpass()

    finder = EquivalentSetFinder("chris", "$dH0sv69")
    finder.find_unique_sets(cols = ["space_time_avg","mass_flow_brush_feeder_avg"])
