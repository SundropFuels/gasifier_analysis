"""LabFunctionLib.py
   General data processing and handling functions for experimental analysis
   Chris Perkins
   2012-01-23

   Version history:
   Version 0.1 - First porting over from R libraries

"""

import numpy as np
import db_toolbox as SQL
import dataFrame_pd as df
import datetime
import unitConversion as uc
import Thermo
import matplotlib.pyplot as plt
import collections
import element_parser as ep
import scipy.optimize as spo
import Cantera as ct

conv = uc.UnitConverter()

class lflException(df.dfException):
    pass

class SQLInterfaceError(lflException):
    pass

class TimeError(lflException):
    pass

class InterpArgumentError(lflException):
    pass

class SQLInterfaceError(lflException):
    pass

class UnitConversionError(lflException):
    pass

class NoUnitsError(UnitConversionError):
    pass

class SpeciesNotDefinedError(lflException):
    pass

class BadStreamError(lflException):
    pass

class ConversionError(lflException):
    pass

class ts_data(df.Dataframe):
    """General timeseries data class"""
    def __init__(self, start, end, data = None, units_dict = None):
        df.Dataframe.__init__(self, data = data, units_dict = units_dict)
        if type(start) !=  datetime.datetime or type(end) !=  datetime.datetime:
            raise TimeError, "The start and end times for a timeseries must be datetime.datetime objects"
        if start > end:
            raise TimeError, "The start time should happen before the end time"

        self.start = start
        self.end = end
        self.uncertainty_dict = {}  #The uncertainty dictionary -- need to add a function to fill this from an SQL table

        self.avgs = {}
        self.stdevs = {}

    def reinitialize_data(self, data):
        self.__init__(self.start, self.end, data = data, units_dict = self.units)

    def SQL_load(self, db_interface, table = ""):
        if not isinstance(db_interface, SQL.db_interface):
            raise SQLInterfaceError, "The passed interface to the database is not valid"

        start_condition = "timestamp >=  '%s'" % self.start.strftime("%Y-%m-%d %H:%M:%S")
        end_condition = "timestamp <=  '%s'" % self.end.strftime("%Y-%m-%d %H:%M:%S")
        try:
            self.SQL_load_data(db_interface, table, conditions = [start_condition, end_condition])
        except SQL.NoDBInUseError:
            raise SQLInterfaceError, "No database selected for the passed interface"
        except SQL.InterfaceNotConnectedError:
            raise SQLInterfaceError, "Passed interface is not connected to a server"
        counter = np.arange(len(self))
        self['counter'] = counter

        for i in self.columns:
            try:
                q=SQL.select_Query(objects=['units'], table='tag_glossary_tbl', condition_list=["simple_name='%s'" % i])
                self.units[i]=db_interface.query(q)[0]['units']
            except IndexError:
                self.units[i]=None

        #Need to add the units load to this and to the unittests

    def SQL_db_upload(self, db_interface, table = ""):
        if not isinstance(db_interface, SQL.db_interface):
            raise SQLInterfaceError, "The passed interface to the database is not valid"

        try:
            self.SQL_upload_data(db_interface, table)
        except SQL.NoDBInUseError:
            raise SQLInterfaceError, "No database selected for the passed interface"
        except SQL.InterfaceNotConnectedError:
            raise SQLInterfaceError, "Passed interface is not connected to a server"
        except df.dfException, e:
            raise SQLInterfaceError, "%s" % e


    def SQL_uc_load(self, db_interface, table = ""):
        pass

    def interpolate_col(self, x, y):
        """Interpolates between None or np.nan values on a given column ('y') relative to 'x' and adds the interpolated column as 'Name'_interp"""
        if x not in self.columns or y not in self.columns:
            raise InterpArgumentError, "The given column name is not in the data frame -- cannot interpolate"

        interp_col = self[y].copy()
        #Scroll forward to find the first non-nan/non-None value
        i = 0
        try:
            while (self[y][i] ==  None or np.isnan(self[y][i])) and i < len(self):
                i+= 1
            if i ==  len(self):
                #The column was empty -- nothing to do
                return interp_col
            begin = i
            while i < len(self):
                i +=  1
                if i ==  len(self):
                    #Trailing nones/nans
                    return interp_col
                elif self[y][i] ==  None or np.isnan(self[y][i]):
                    continue
                elif i ==  len(self):
                    #Trailing nones/nans
                    return interp_col
                elif i - begin ==  1:
                    #There were no nans or Nones in the mix, move forward
                    begin = i
                else:
                    #Now interpolate between the values
                    end = i
                    for j in range(begin+1,end):
                        interp_col[j] = (self[y][end] - self[y][begin])/(self[x][end] - self[x][begin])*(self[x][j] - self[x][begin]) + self[y][begin]
                    begin = i
        except IndexError:
            return interp_col
        return interp_col
                        
    
    def get_ms_from_csv(self, csvfile):
        ppmlist = ['Benzene_MS', 'Hydrogen Sulfide_MS', 'Napthalene_MS', 'Toluene_MS']
        for i in self.data:
            if i.endswith('_MS'):
                for j in range(len(self[i])): #Clear list
                    self[i][j] = None
                with open(csvfile) as f:
                    colnames = f.readline().split(',')
                    timecol = colnames.index('Time&Date')
                    icol = colnames.index(i.replace('_MS', ''))
                    for line in f:
                        row = line.split(',')
                        time = datetime.datetime.strptime(row[timecol],'%Y/%m/%d %H:%M:%S')
                        if i in ppmlist:
                            conc = float(row[icol])/10000
                        else:
                            conc = float(row[icol])
                        self[i][np.where(self['timestamp'] == time)[0]] = conc
                    
    def calc_rel_time(self):
        rel_time = []
        for i in self['timestamp']:
            rel_time.append((i-self['timestamp'][0]).seconds)
        rel_time = np.array(rel_time)
        return rel_time

    def rid_nones(self, colname):
        for i in range(len(colname)-1):
            if colname[i-1] == None and colname[i]!= None:
                firstval = i
            if colname[i+1] == None and colname[i]!= None:
                lastval = i
        colname[0:firstval] = colname[firstval]
        colname[lastval+1:] = colname[lastval]
        return(colname)

    def generate_averages_and_stdevs(self, cols=None):
        """Generates numbers for the averages and standard deviations of the indicated columns, and stores them in member dicts"""
        if cols is None:
            cols = self.data.keys()

        for key in cols:
            try:
                #by default, I ignore nan and inf values
                
                self.avgs[key] = self[key][np.isfinite(self[key].astype(float))].mean()
                self.stdevs[key] = self[key][np.isfinite(self[key].astype(float))].std()
            except KeyError:
                raise lflExeception, "%s is not a key in the dataframe" % key
            except ZeroDivisionError:
                self.avgs[key] = np.nan
                self.stdevs[key] = np.nan
	
    def interp_ga(self):
        self['rel_time'] = self.calc_rel_time()
        self.units['rel_time']='s'
        for i in self.data:
            if i.endswith('_MS') or i.endswith('_GC'):
                colname = i
                try:
                    interpcol = self.interpolate_col('rel_time', colname)
                    interpcol = self.rid_nones(interpcol)
                    self[colname] = interpcol
                except IndexError: pass
                if self.units[i]=='ppm':
                    self[i]=self[i]/10000
                    self.units[i]='%'
                
   
    def calc_product_flowrates(self):
        """Calculate product gas flow rates in moles/s"""
        gaslist = []
        for i in self.data:
            if i.endswith('_MS'):
                gaslist.append(i)
        for i in gaslist:
            product_flowrate = self.data[i]/100*self.data['outlet_flowrate']
            self.data[i+'_flowrate'] = product_flowrate

    
    def calc_inst_conversion(self):
        """Calculates instantaneous Carbon conversion"""
        self['carbon_conversion'] = self['carbon_out_total']/self['carbon_in']

    

    def calc_average_feed_rate(self):
        return np.average(self.data['mass_flow_brush_feeder'])

    def calc_std_feed_rate(self):
        return np.std(self.data['mass_flow_brush_feeder'])
        
    def calc_total_c_in(self):
        return np.sum(self.data['carbon_in'])

    def calc_total_c_out(self):
        return np.sum(self.data['carbon_out_total'])

    def calc_carbon_x(self):
        return np.sum(self.data['carbon_out_total'])/np.sum(self.data['carbon_in'])

    

    def convert_col_units(self, column, units):
        """Will convert the given column into the new units from the existing units"""
        
        if column not in self.columns:
            raise df.NoColumnError, "The requested column is not in the data frame."

        if self.units ==  {}:
            raise NoUnitsError, "The units have not been defined for this column"

        conv = uc.UnitConverter()

        try:
            self[column] = conv.convert_units(self[column], self.units[column], units)
        except uc.UnitNotFoundError:
            raise UnitConversionError, "The unit was not in the database"
        except uc.InconsistentUnitError:
            raise UnitConversionError, "The units are not dimensionally consistent with the current units"

        self.units[column] = units    

    def get_val(self, col_name):
        if col_name in self.units.keys():
            return (self[col_name], self.units[col_name])
        else:
            return (self[col_name], "")

    def list_keys(self):
        """Generates list of column names in their respective categories"""
        keylist = collections.OrderedDict({})
        keylist['MFCs'] = []
        keylist['Temperatures'] = []
        keylist['Pressures'] = []
        keylist['Setpoints'] = []
        keylist['Mass Spec'] = []
        keylist['GC'] = []
        keylist['NDIR'] = []
        keylist['Product Flow Rates'] = []
        keylist['Carbon Flows'] = []
        
        for i in self.data:
            if i.endswith('_MS'):
                keylist['Mass Spec'].append(i)
            if i.endswith('_GC'):
                keylist['GC'].append(i)
            if i.startswith('temp_'):
                keylist['Temperatures'].append(i)
            if i.startswith('pressure_'):
                keylist['Pressures'].append(i)
            if i.startswith('mass_'):
                keylist['MFCs'].append(i)
            if i.endswith('_flowrate'):
                keylist['Product Flow Rates'].append(i)
            if i.startswith('carbon_') or i.endswith('_carbon_out'):
                keylist['Carbon Flows'].append(i)
            if i.startswith('setpoint_'):
                keylist['Setpoints'].append(i)
            if i.startswith('NDIR_'):
                keylist['NDIR'].append(i)
        for k in keylist:
            keylist[k].sort()
        return keylist

class Biomass:
    def __init__(self, sample = None, CHN = None, moisture = None, enthalpy_formation = None):
        #Heat of formation SEP SYP -3.156 MJ/kg or -1432 kJ/lb
        self.sample=sample
        if CHN is None:
            CHN = {}
        self.CHN = CHN #{C:(Value, 'Units'), H:(Value, 'Units'), N:(Value, 'Units')}
        self.moisture = moisture #(Value, 'Units')
        self.enthalpy_formation = enthalpy_formation #(Value, 'Units')

    def lookup(self):
        pass

    def calc_CHN_dry(self):
        self.CHN_dry={}
        self.CHN_dry['C']=(conv.convert_units(self.CHN['C'][0], self.CHN['C'][1], 'fraction')/
                           (1.-conv.convert_units(self.moisture[0], self.moisture[1], 'fraction')), 'fraction')
        self.CHN_dry['H']=((conv.convert_units(self.CHN['H'][0], self.CHN['H'][1], 'fraction')-
                            conv.convert_units(self.moisture[0], self.moisture[1], 'fraction')*
                           2.02/18.02)/(1.-conv.convert_units(self.moisture[0], self.moisture[1], 'fraction')), 'fraction')
        self.CHN_dry['N']=(conv.convert_units(self.CHN['N'][0], self.CHN['N'][1], 'fraction')/
                           (1.-conv.convert_units(self.moisture[0], self.moisture[1], 'fraction')), 'fraction')

class Gasification_Experiment:
    def __init__(self, name, biomass_sample = None, biomass_feedrate = None,
                 entrainment_flow = None, entrainment_gas = None,
                 makeup_flow = None, makeup_gas = None,
                 argon_flow = None, fill_flow = None, temperature = None,
                 pressure = None, start_time = None, end_time = None,
                 steam_flow = None):
        self.biomass_sample = biomass_sample
        self.biomass_feedrate = biomass_feedrate #(value, units)
        self.entrainment_flow = entrainment_flow #(value, units)
        self.entrainment_gas = entrainment_gas
        self.makeup_flow = makeup_flow #(value, units)
        self.makeup_gas = makeup_gas
        self.argon_flow = argon_flow #(value, units)
        self.fill_flow = fill_flow #(value, units)
        self.temperature = temperature #(value, units)
        self.pressure = pressure #(value, units)
        self.steam_flow = steam_flow #(value, units)
        self.start_time = start_time
        self.end_time = end_time
        
    def lookup(self):
        pass
        
class Stream:
    species = {}
    species['CO'] = {'C':1,'O':1, 'name':'Carbon Monoxide'}
    species['CO2'] = {'C':1,'O':2, 'name':'Carbon Dioxide'}
    species['CH4'] = {'C':1, 'H':4, 'name':'Methane'}
    species['H2'] = {'H':2, 'name':'Hydrogen'}
    species['H2S'] = {'H':2, 'S':1, 'name':'Hydrogen Sulfide'}
    species['C2H2'] = {'C':2,'H':2, 'name':'Acetylene'}
    species['C2H4'] = {'C':2,'H':4, 'name':'Ethylene'}
    species['C2H6'] = {'C':2,'H':6, 'name':'Ethane'}
    species['C3H8'] = {'C':3, 'H':8, 'name':'Propane'}
    species['C3H6'] = {'C':3, 'H':6, 'name':'Propylene'}
    species['C6H6'] = {'C':6, 'H':6, 'name':'Benzene'}
    species['C7H8'] = {'C':7, 'H':8, 'name':'Toluene'}
    species['C10H8'] = {'C':10, 'H':8, 'name':'Napthalene'}
    species['H2O'] = {'H':2,'O':1, 'name':'Water'}
    species['Ar'] = {'Ar':1, 'name':'Argon'}
    species['N2'] = {'N':2, 'name':'Nitrogen'}
    species['CELL'] = {'C':6, 'H':10, 'O':5, 'name':'Cellulose'}
    species['HCE'] = {'C':5, 'H':8, 'O':4, 'name':'Hemicellulose'}
    species['LIGH'] = {'C':22, 'H':28, 'O':9, 'name':'Lig-H'}
    species['LIGO'] = {'C':20, 'H':22, 'O':10, 'name':'Lig-O'}
    species['LIGC'] = {'C':15, 'H':14, 'O':4, 'name':'Lig-C'}
    
    ct_trans = {}
    ct_trans['Ar'] = {'AR':1.0}
    ct_trans['biomass'] = {'CELL':0.3932, 'LIGC':0.10, 'LIGH':0.10, 'LIGO':0.10, 'HCE':0.30}   #Still missing ASH - will correct in a bit

    names = {}
    for i in species.keys():       
        names[species[i]['name']] = i

    def __init__(self, name, flowrate = None, composition = None,
                 basis = "molar", temperature = None, pressure = None,
                 density = None, compressible = None, std_temperature = (25.0, 'C'), std_pressure = (101325.0, 'Pa')):

        if composition is None:
            composition = {}
        self.name = name
        self.composition = composition #{specie:fraction}
        self.basis = basis
        self.temperature = temperature #(value, units)
        self.pressure = pressure #(value, units)
        self.density = density #(value, units)
        self.flowrate = flowrate #(value, units)
        self.compressible = compressible
        self.std_temperature = std_temperature
        self.std_pressure = std_pressure

        self.special_species = {}
        
        #Create Cantera object for the stream
        self.ctphase = ct.importPhase('cantera_biomass/GasifierSpecies.cti','gas')
        if self.temperature is not None:
            self.ctphase.set(T = conv.convert_units(self.temperature[0], self.temperature[1], 'K'))
        if self.pressure is not None:
            self.ctphase.set(P = conv.convert_units(self.pressure[0], self.pressure[1], 'kg/s^2/m'))
        if self.composition is not None:
            self.ct_setcomp()
        #will need to set up None checking in the enthalpy function to make sure that stream values are not empty    
        
    def gas_volumetric_flowrate(self, units):
        """Returns the gas volumetric flowrate, in the desired units"""
        conv = uc.UnitConverter()
        if self.basis == "gas_volume":
            return conv.convert_units(self.flowrate[0], self.flowrate[1], units)
        elif self.basis == "std_gas_volume":
            T = conv.convert_units(self.temperature[0], self.temperature[1], 'K')
            p = conv.convert_units(self.pressure[0], self.pressure[1], 'Pa')
            std_T = conv.convert_units(self.std_temperature[0], self.std_temperature[1], 'K')
            std_P = conv.convert_units(self.std_pressure[0], self.std_pressure[1], 'Pa')
            return conv.convert_units(self.flowrate[0], self.flowrate[1], units)*T/std_T*std_P/p
        elif self.basis == "molar":
            f = conv.convert_units(self.flowrate[0], self.flowrate[1], 'mol/s')
            T = conv.convert_units(self.temperature[0], self.temperature[1], 'K')
            p = conv.convert_units(self.pressure[0], self.pressure[1], 'Pa')
            return conv.convert_units(8.314*f*T/p, 'm^3/s', units)
        elif self.basis == "mass":
            #convert to a molar flow first ###!!!###
            val = conv.convert_units(self.flowrate[0], self.flowrate[1], 'g/s')
            for species in self.composition.keys():
                if species in SpecialMolecule.MW.keys():
                    MW = SpecialMolecule.MW[species]
                else:
                    MW = 0
                    breakdown = ep.parse_species(species)
                    try:
                        for ele, v in breakdown.items():
                            MW +=  v*Element.MW[ele]
                    except KeyError:
                        raise BadStreamError, "%s does not have an entry in the Element molecular weight dictionary" % ele
            return val*self.composition[species]/MW


    def elementalFactor(self, element):
        """Returns the units of element/basis unit of a given element in the feed"""
        factor = 0.0
        for (specie, fraction) in self.composition.items():
            try:
                spe_dict = Stream.species[specie]
                
            except KeyError:
                try:
                    spe_dict = self.special_species[specie]
                    
                except KeyError:
                    raise SpeciesNotDefinedError, "%s does not have an elemental breakdown definition...check special_species{}" % specie
            try:
                
                factor +=  fraction*spe_dict[element]
            except KeyError:
                pass #It's okay to pass here, because not all species will have the specified element
        
        return factor

    def speciesFactor(self, specie):
        """Returns the units of species/basis unit of a given element in the feed"""
        try:
            return self.composition[specie]
        except KeyError:
            pass

    def calcElementalMolarFlowrate(self, element):
        """Calculates the molar feed rates (mol/s) of an element for a given stream flowrate"""
        if self.basis ==  "molar":
            #First need to convert to mol/s
            conv = uc.UnitConverter()
            val = conv.convert_units(self.flowrate[0], self.flowrate[1], 'mol/s')
            return val*self.elementalFactor(element)
        elif self.basis ==  "mass":
            #First need to convert to g/s
            conv = uc.UnitConverter()
            val = conv.convert_units(self.flowrate[0], self.flowrate[1], 'g/s')
            
            return val*self.elementalFactor(element)/Element.MW[element]
        elif self.basis ==  "gas_volume":
            #First need to convert to m^3/s
            conv = uc.UnitConverter()
            val = conv.convert_units(self.flowrate[0], self.flowrate[1], 'm^3/s')
            p = conv.convert_units(self.pressure[0], self.pressure[1], 'kg/s^2/m')
            T = conv.convert_units(self.temperature[0], self.temperature[1], 'K')
            #Need to convert the temperature to K
            return val*p/(8.314*T)*self.elementalFactor(element)
        elif self.basis == "std_gas_volume":
            #First need to convert to Sm^3/s
            conv = uc.UnitConverter()
            val = conv.convert_units(self.flowrate[0], self.flowrate[1], 'm^3/s')
            p = conv.convert_units(self.std_pressure[0], self.std_pressure[1], 'kg/s^2/m')
            T = conv.convert_units(self.std_temperature[0], self.std_temperature[1], 'K')
            return val*p/(8.314*T)*self.elementalFactor(element)

        elif self.basis ==  "liquid_volume":
            #First need to convert the density to the same units as the liquid
            pass
        else:
            pass

    def calcSpeciesMolarFlowrate(self, species):
        """Calculates the molar feed rate (mol/s) of a species for a given stream flowrate"""
        if species not in self.composition.keys():
            return np.zeros(len(self.flowrate[0]))

        if self.basis ==  "molar":
            conv = uc.UnitConverter()
            val = conv.convert_units(self.flowrate[0], self.flowrate[1], 'mol/s')
            return val*self.composition[species]
        elif self.basis ==  "mass":
            #Need to convert to molar amounts based on the molecular weight
            conv = uc.UnitConverter()
            val = conv.convert_units(self.flowrate[0], self.flowrate[1], 'g/s')
            
            if species in SpecialMolecule.MW.keys():
                MW = SpecialMolecule.MW[species]
            else:
                MW = 0
                breakdown = ep.parse_species(species)
                try:
                    for ele, v in breakdown.items():
                        MW +=  v*Element.MW[ele]
                except KeyError:
                    raise BadStreamError, "%s does not have an entry in the Element molecular weight dictionary" % ele
            return val*self.composition[species]/MW
        elif self.basis ==  "gas_volume":
            
            conv = uc.UnitConverter()
            val = conv.convert_units(self.flowrate[0], self.flowrate[1], 'm^3/s')
            p = conv.convert_units(self.pressure[0], self.pressure[1], 'kg/s^2/m')
            T = conv.convert_units(self.temperature[0], self.temperature[1], 'K')
            return val*self.composition[species]*p/(8.314*T)
        elif self.basis == "std_gas_volume":
            conv = uc.UnitConverter()
            val = conv.convert_units(self.flowrate[0], self.flowrate[1], 'm^3/s')
            
            p = conv.convert_units(self.std_pressure[0], self.std_pressure[1], 'Pa')
            T = conv.convert_units(self.std_temperature[0], self.std_temperature[1], 'K')
            return val*self.composition[species]*p/(8.314*T)

        else:
            raise lflException, "Basis not recognized!"


    def calcSpeciesVolumetricFlowrate(self, species):    
        """Calculates the volumetric feed rate (m^3/s) of a species for a given stream flowrate"""
        #This needs to be generalized to admit liquid states -- also, unit tests must be written!
        if species not in self.composition.keys():
            return np.zeros(len(self.flowrate[0]))

        if self.basis ==  "molar":
            conv = uc.UnitConverter()
            val = conv.convert_units(self.flowrate[0], self.flowrate[1], 'mol/s')
            T = conv.convert_units(self.temperature[0], self.temperature[1], 'K')
            P = conv.convert_units(self.pressure[0], self.pressure[1], 'kg/s^2/m')
            return val*self.composition[species]*8.314*T/P
        elif self.basis ==  "mass":
            #Need to convert to molar amounts based on the molecular weight
            conv = uc.UnitConverter()
            val = conv.convert_units(self.flowrate[0], self.flowrate[1], 'g/s')
            
            if species in SpecialMolecule.MW.keys():
                MW = SpecialMolecule.MW[species]
            else:
                MW = 0
                breakdown = ep.parse_species(species)
                try:
                    for ele, v in breakdown.items():
                        MW +=  v*Element.MW[ele]
                except KeyError:
                    raise BadStreamError, "%s does not have an entry in the Element molecular weight dictionary" % ele
            molar = val*self.composition[species]/MW
            T = conv.convert_units(self.temperature[0], self.temperature[1], 'K')
            P = conv.convert_units(self.pressure[0], self.pressure[1], 'kg/s^2/m')
            return molar * T*8.314/P
        elif self.basis ==  "gas_volume":
            conv = uc.UnitConverter()
            val = conv.convert_units(self.flowrate[0], self.flowrate[1], 'm^3/s')          
            return val*self.composition[species]
        elif self.basis == "std_gas_volume":
            conv = uc.UnitConverter()
            val = conv.convert_units(self.flowrate[0], self.flowrate[1], 'm^3/s')
            T = conv.convert_units(self.temperature[0], self.temperature[1], 'K')
            p = conv.convert_units(self.pressure[0], self.pressure[1], 'Pa')
            T_std = conv.convert_units(self.std_temperature[0], self.std_temperature[1], 'K')
            p_std = conv.convert_units(self.std_pressure[0], self.std_pressure[1], 'Pa')
            val *= T/T_std * p_std/p
            return val*self.composition[species]

        else:
            raise lflException, "Basis not recognized!"

    def _parse_species(self, species_str):
        """Parses a molecule into into a dictionary of {element: # of atoms}"""

        parsed_elements = {}
        current_element = ""
        current_number_str = ""
    

        i = 0
        if not species_str[0].isalpha():
            raise BadCharacterError, "A molecule must start with an alphabetical character"

        while i < len(species_str):
            current_char = species_str[i]
        
            if current_char.isalpha():
                if i+1 ==  len(species_str) or species_str[i+1].isupper():
                    #Allows for single character names, like CH4
                    current_element = "".join([current_element, current_char])
                    if current_element in parsed_elements.keys():
                        parsed_elements[current_element] +=  1
                    else:
                        parsed_elements[current_element] = 1
                    current_element = ""
                
                
                else:
                
                    current_element = "".join([current_element, current_char])
                i +=  1
                continue

            elif current_char.isdigit():
                #we have gotten to the end of an element name
                print current_element
                while i < len(species_str) and species_str[i].isdigit():
                    print current_char
                    current_char = species_str[i]
                    current_number_str = "".join([current_number_str, current_char])
                    i +=  1
                if current_number_str ==  '':
                    raise BadCharacterError, "Each element must have a number associated with it"
                if current_element in parsed_elements.keys():
                    parsed_elements[current_element] +=  int(current_number_str)
                else:
                    parsed_elements[current_element] = int(current_number_str)
            
                current_element = ""
                current_number_str = ""
            else:
                raise BadCharacterError, "A molecule can only contain alphabetical and numerical characters"

        return parsed_elements

    def get_enthalpy(self, units):
        if not IsInstance(units, str):
            raise UnitConversionError, "The provided units must be a string for the enthalpy function"
        
        self._calc_enthalpy()
        conv = uc.UnitConverter()
        return conv.convert_units(self.enthalpy[0], self.enthalpy[1], units)
       

    def get_entropy(self, units):
        if not IsInstance(units, str):
            raise UnitConversionError, "The provided units must be a string for the enthalpy function"
        
        self._calc_entropy()
        conv = uc.UnitConverter()
        return conv.convert_units(self.entropy[0], self.entropy[1], units)

    def ct_setcomp(self):
        """Returns a string that can be used to input the stream's composition into a Cantera object"""
        #First, need to format everything into the Cantera Species
        if self.basis == 'molar' or self.basis == 'gas_volume' or self.basis == 'std_gas_volume':
            #This is the easiest case -- just pull all the compositional values and set the phase appropriately
            set_string = ""
            for specie in self.composition.keys():
                if specie in Stream.ct_trans:
                    for sub in Stream.ct_trans[specie]:
                        set_string += '%s:%s' % (sub, Stream.ct_trans[specie][sub]*self.composition[specie])  #This, of course, assumes that the basis in ct_trans is the same as in self.composition
                else:
                    set_string += '%s:%s' % (specie, self.composition[specie])
            #set the phase composition
            self.ctphase.setMoleFractions(set_string)

        elif self.basis == 'mass':
            #This is also pretty easy -- just pull all the compositional values and set the phase appropriately
            set_string = ""
            for specie in self.composition.keys():
                if specie in Stream.ct_trans:
                    for sub in Stream.ct_trans[specie]:
                        set_string += '%s:%s' % (sub, Stream.ct_trans[specie][sub]*self.composition[specie])  #This, of course, assumes that the basis in ct_trans is the same as in self.composition
                else:
                    set_string += '%s:%s' % (specie, self.composition[specie])
            #set the phase composition
            self.ctphase.setMassFractions(set_string)

        else:
            raise lflException, '%s is not a valid stream basis' % self.basis
        
        
    def _calc_enthalpy(self):
        """Calculates the stream enthalpy and stores it in self.enthalpy"""
        
        if self.temperature==None:
            raise BadStreamError, 'Stream temperature is not defined.'
        if self.pressure==None:
            raise BadStreamError, 'Stream pressure is not defined.'
        conv = uc.UnitConverter()
        if self.basis == 'molar':
            #convert to kmol/s:
            flow = conv.convert_units(self.flowrate[0], self.flowrate[1], 'kmol/s')
            self.enthalpy = [flow*self.ctphase.enthalpy_mole(), 'J/s']

        elif self.basis == 'mass':
            #convert to kg/s
            flow = conv.convert_units(self.flowrate[0], self.flowrate[1], 'kg/s')
            self.enthalpy = [flow*self.ctphase.enthalpy_mass(), 'J/s']

        elif self.basis ==  "gas_volume":
                        
            val = conv.convert_units(self.flowrate[0], self.flowrate[1], 'm^3/s')
            p = conv.convert_units(self.pressure[0], self.pressure[1], 'kg/s^2/m')
            T = conv.convert_units(self.temperature[0], self.temperature[1], 'K')
            flow =  val*p/(8.314*T)
            self.enthalpy = [flow*self.ctphase.enthalpy_mole(), 'J/s']

        elif self.basis == "std_gas_volume":
            
            val = conv.convert_units(self.flowrate[0], self.flowrate[1], 'm^3/s')
            
            p = conv.convert_units(self.std_pressure[0], self.std_pressure[1], 'Pa')
            T = conv.convert_units(self.std_temperature[0], self.std_temperature[1], 'K')
            flow =  val*p/(8.314*T)
            self.enthalpy = [flow*self.ctphase.enthalpy_mole(), 'J/s']
            

    def _calc_entropy(self):
        """Calculates the stream entropy and stores it in self.entropy"""
        conv = uc.UnitConverter()
        if self.temperature==None:
            raise BadStreamError, 'Stream temperature is not defined.'
        if self.pressure==None:
            raise BadStreamError, 'Stream pressure is not defined.'
        if self.basis == 'molar':
            #convert to kmol/s:
            flow = conv.convert_units(self.flowrate[0], self.flowrate[1], 'kmol/s')
            self.entropy = [flow*self.ctphase.entropy_mole(), 'J/K/s']

        elif self.basis == 'mass':
            #convert to kg/s
            flow = conv.convert_units(self.flowrate[0], self.flowrate[1], 'kg/s')
            self.entropy = [flow*self.ctphase.entropy_mass(), 'J/K/s']

        elif self.basis ==  "gas_volume":
                        
            val = conv.convert_units(self.flowrate[0], self.flowrate[1], 'm^3/s')
            p = conv.convert_units(self.pressure[0], self.pressure[1], 'kg/s^2/m')
            T = conv.convert_units(self.temperature[0], self.temperature[1], 'K')
            flow =  val*p/(8.314*T)
            self.entropy = [flow*self.ctphase.entropy_mole(), 'J/K/s']

        elif self.basis == "std_gas_volume":
            
            val = conv.convert_units(self.flowrate[0], self.flowrate[1], 'm^3/s')
            
            p = conv.convert_units(self.std_pressure[0], self.std_pressure[1], 'Pa')
            T = conv.convert_units(self.std_temperature[0], self.std_temperature[1], 'K')
            flow =  val*p/(8.314*T)
            self.entropy = [flow*self.ctphase.entropy_mole(), 'J/K/s']
        
    def convert_to_mass_basis(self):
        if self.basis == 'mass':
            #do nothing
            pass
        
        elif self.basis == 'gas_volume':
            #I'm going to be sneaky here and convert to molar first, then let the execution drop downward
            conv = uc.UnitConverter()
            p = conv.convert_units(self.pressure[0], self.pressure[1], 'Pa')
            T = conv.convert_units(self.temperature[0], self.temperature[1], 'K')
            self.flowrate = [conv.convert_units(self.flowrate[0], self.flowrate[1], 'm^3/s')*p/(8.314*T), 'mol/s']
            self.basis = 'molar'

        elif self.basis == 'std_gas_volume':
            conv = uc.UnitConverter()
            p = conv.convert_units(self.std_pressure[0], self.std_pressure[1], 'Pa')
            T = conv.convert_units(self.std_temperature[0], self.std_temperature[1], 'K')
            self.flowrate = [conv.convert_units(self.flowrate[0], self.flowrate[1], 'm^3/s')*p/(8.314*T), 'mol/s']
            self.basis = 'molar'

        if self.basis == 'molar':
            #calculate the molar compositions from the mass compositions
            
            conv = uc.UnitConverter()
            #calculate the average molecular weight first
            species_dict = {}
            for species in self.composition:
            
                if species in SpecialMolecule.MW.keys():
                    MW = SpecialMolecule.MW[species]
                else:
                    MW = 0
                    breakdown = ep.parse_species(species)
                    try:
                        for ele, v in breakdown.items():
                            MW +=  v*Element.MW[ele]
                    except KeyError:
                        raise BadStreamError, "%s does not have an entry in the Element molecular weight dictionary" % ele
                species_dict[species] = MW*self.composition[species]
            avg_MW = sum(species_dict.values())
            for species in species_dict:
                self.composition[species] = species_dict[species]/avg_MW
            self.flowrate = [avg_MW*conv.convert_units(self.flowrate[0],self.flowrate[1], 'mol/s'), 'g/s']
            self.basis = 'mass'

        
            



class ProcessObject:
    
    def __init__(self, inlets, outlets = None):
        
        #Inlet Checks.  Maybe in future don't need to have inlets fully defined if outlets are, but not sure how to implement now.
        if type(inlets) != list:
            raise BadStreamError, 'Inlets must be input as a list of Stream objects.'
        for inlet in inlets:
            if not isinstance(inlet, Stream):
                raise BadStreamError, 'Inlet %s is not a Stream object.' %inlet.name
            #I'm going to relax this condition -- maybe we want the ProcessObject to SOLVE for T, P, or other
            """
            if inlet.temperature = None:
                raise BadStreamError, 'Inlet %s temperature is not defined.' %inlet.name
            if inlet.pressure = None:
                raise BadStreamError, 'Inlet %s pressure is not defined.' %inlet.name
            if inlet.composition = None:
                raise BadStreamError, 'Inlet %s composition is not defined.' %inlet.name
            if inlet.flowrate = None:
                raise BadStreamError, 'Inlet %s flow rate is not defined.' %inlet.name
            """
            
                
        #Outlet Checks
        if outlets is not None:
            if type(outlets) != list:
                raise BadStreamError, 'Outlets must be input as a list of Stream objects.'
            for outlet in outlets:
                if not isinstance(outlet, Stream):
                    raise BadStreamError, 'Outlet %s is not a Stream object.' %outlet.name

        #Probably don't want this as general behavior, either
        """        
        #Find total flowrates for all species in inlets
        self.inspeciesflowrates={}        
        for inlet in self.inlets:
            for sp in inlet.composition.keys():
                if sp not in self.inspeciesflowrates:
                    self.inspeciesflowrates[sp]=0
                self.inspeciesflowrates[sp]+=inlet.calcSpeciesMolarFlowrate(sp)
       
        """
    #There has GOT to be a way to generalize this with decorators or similar
    def totalInletEnthalpy(self, units):
        H = 0
        for stream in self.inlets:
            H += stream.get_enthalpy(units)
        return H

    def totalOutletEnthalpy(self, units):
        H = 0
        for stream in self.outlets:
            H += stream.get_enthalpy(units)
        return H

    def totalInletEntropy(self,units):
        S = 0
        for stream in self.inlets:
            S += stream.get_entropy(units)
        return S

    def totalOutletEntropy(self,units):
        S = 0
        for stream in self.outlets:
            S += stream.get_entropy(units)
        return S

    def deltaH(self, units):
        return self.totalOutletEnthalpy(units) - self.totalInletEnthalpy(units)
    
    def deltaS(self, units):
        return self.totalOutletEntropy(units) - self.totalInletEnthalpy(units)
            
class Mixer(ProcessObject):
    
    """A mixer blends multiple inlet streams into one outlet stream"""

    def __init__(self, name, outlet_pressure = None, **kwargs):
        ProcessObject.__init__(**kwargs)
        
        self.outlets = [Stream(name = '%s_outlet' % name)]
        
        #Need to solve for the outlet stream pressure, temperature, and composition
        #pressure is easy -- it is assumed that the pressure drops to the LOWEST stream pressure entering the mixer if an outlet pressure is not specified
        self._calc_outlet_pressure(outlet_pressure)
        self._calc_outlet_flowrate()
        self._calc_outlet_temperature()

    def recalc(self, outlet_pressure = None):
        self._calc_outlet_pressure(outlet_pressure)
        self._calc_outlet_flowrate()
        self._calc_outlet_temperature()

    def _calc_outlet_pressure(self, outlet_pressure):
        conv = uc.UnitConverter()
        minP = None
        if outlet_pressure is None:
            for stream in self.inlets:
                if minP is None or conv.convert_units(stream.pressure[0], stream.pressure[1], 'Pa') < minP:
                    minP = conv.convert_units(stream.pressure[0], stream.pressure[1], 'Pa')
            else:
                minP = outlet_pressure
        self.outlets[0].pressure = [minP, 'Pa']

    def _calc_outlet_flowrate(self):
        conv = uc.UnitConverter()
        #Need to put everything on a consistent basis - if everything the same, just add them all together; otherwise use mass
        #Also, need to set the composition of the outlet
        c = True
        basis_fl_dict = {'molar':'mol/s', 'mass':'kg/s', 'gas_volume':'m^3/s', 'std_gas_volume':'m^3/s'}
        basis_choice = self.inlets[0].basis
        for inlet in self.inlets:
            if inlet.basis != basis_choice or inlet.basis == 'gas_volume' or inlet.basis == 'std_gas_volume':
                c = False

        if not c:
            #convert all streams to a mass basis
            for inlet in self.inlets:
                if inlet.basis != 'mass':
                    inlet.convert_to_mass_basis()
            basis_choice = 'mass'

        
        self.outlets[0].basis = basis_choice
        fl_sum = 0
        for inlet in self.inlets:
            fl_sum += conv.convert_units(inlet.flowrate[0], inlet.flowrate[1], basis_fl_dict[basis_choice])
        self.outlets[0].flowrate = (fl_sum, basis_fl_dict[basis_choice])
        #need to generate a total species list for compositional matching - then composition is simply the sum of the streams for that species divided by the total flowrate
        species_list = []
        for inlet in self.inlets:
            for species in inlet.composition: 
                if species not in species_list:
                    species_list.append(species)
        composition = {}
        for specie in species_list:
            spec_sum = 0
            for inlet in self.inlets:
                spec_sum += conv.convert_units(inlet.flowrate[0], inlet.flowrate[1], basis_fl_dict[basis_choice])*inlet.composition[species]
            composition[species] = spec_sum/fl_sum
        self.outlets[0].composition = composition

        #Just driving all volume directly to mass now -- should work


    def _calc_outlet_temperature(self):
        #Solving the equation dH = 0 (not a heat exchanger, so Q and W are both 0)
        #guess a temperature for the outlet as a mean of the inlet temperatures
        conv = uc.UnitConverter()
        temp_sum = 0.0
        for inlet in self.inlets:
            temp_sum += conv.convert_units(inlet.temperature[0], inlet.temperature[1], 'K')
        temp_avg = temp_sum/len(self.inlets)
        outlet_temp = spo.newton(func = self.deltaH, x0 = temp_avg, args = ('J/s'))
        self.outlets[0].temperature = [outlet_temp, 'K'] 



    
        
class Reactor(ProcessObject):
    """Reactor Class..."""
    def __init__(self, temperature = None, pressure = None, **kwargs):
        ProcessObject.__init__(**kwargs)
    def calc_species_generation(self):
        pass
    def calc_species_consumption(self):
        pass
    def calc_enthalpy_change(self, units):
        return self.deltaH(units)
    def calc_entropy_change(self, units):
        return self.deltaS(units)
        
class Condenser(ProcessObject):
    def __init__(self, **kwargs):
        ProcessObject.__init__(**kwargs)
        
class PhysicalConstants:
    """Class to hold global physical constant variables"""
    std_temp = 298   #In K
    std_pressure = 100000 #In Pa

class Element:
    """A physical element.  Will have properties of the elements, plus a handy reference library"""
    Elements = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au','Hg', 'Tl', 'Pb','Bi', 'Po', 'At', 'Rn', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu']
    MW = {}
    MW['C'] = 12.01
    MW['H'] = 1.01
    MW['O'] = 16.00
    MW['N'] = 14.01
    MW['S'] = 32.065
    MW['Cl'] = 35.453
    MW['Ar'] = 39.948


class Molecule:
    def __init__(self, formula):
        self.formula = formula
        self.form_dict = ep.parse_species(self.formula)


    def MW(self):
        MW = 0.0
        
        for element in self.form_dict:
            MW += Element.MW[element] * self.form_dict[element]
        return MW

class SpecialMolecule(Molecule):
    """This class holds specially designated molcules, so that I can calculate molecular weights without too much of a problem"""
    MW = {}
    MW['LIGH'] = 436
    MW['LIGO'] = 422
    MW['LIGC'] = 258
    MW['CELL'] = 162
    MW['HCE'] = 132

class ProcTS(ts_data):

    """Timeseries data for chemical processes"""
    def __init__(self, start, end, data = None, units_dict = None):
        ts_data.__init__(self, start, end, data, units_dict)
        
        self.inlet_streams = []
        self.outlet_streams = []
        self.proc_elements = []
        self.proc_species = []
        self.inert_species = []

    

    def generate_inlet_outlet_elemental_flows(self, name_qualifier = None):
        """Generates the standard inlet and outlet elemental flows for elements in self.proc_elements"""
        if self.inlet_streams ==  [] or self.outlet_streams ==  []:
            raise BadStreamError, "Inlet/outlet streams not yet defined"
        
        if name_qualifier ==  None:
            outlet_name = 'outlet'
            inlet_name = 'inlet'

        else:
            if type(name_qualifier) !=  str:
                raise lflException, "The new column name qualifier must be a string"
            else:
                outlet_name = 'outlet_%s' % name_qualifier
          

        
        for element in self.proc_elements:
            if element not in Element.Elements:
                raise BadStreamError, "The given element: %s is not a physical element" % element

            self['%s_%s' % (element,inlet_name)] = self.calc_elemental_molar_feedrates(self.inlet_streams, element)
            self['%s_%s' % (element,outlet_name)] = self.calc_elemental_molar_feedrates(self.outlet_streams, element)
            self.units['%s_%s' % (element,inlet_name)] = 'mol/s'
            self.units['%s_%s' % (element,outlet_name)] = 'mol/s'


    def calc_elemental_molar_feedrates(self, stream_list, element):
        """Returns the elemental molar feedrate in a given set of streams"""
        tot = 0        
        

        for stream in stream_list:
            
            if not isinstance(stream, Stream):
                
                raise BadStreamError, "%s is not a stream, which is required for elemental calculations" % stream
            
            tot +=  stream.calcElementalMolarFlowrate(element)

        return tot

    def outlet_stream_from_tracer(self, inlet_streams, basis, tracer_species, outlet_tracer_concentration, name='outlet'):
        """Generates an outlet stream, on a determined basis, from a mass balance over a list of inlet streams on the tracer species"""
        #inlet_streams is a list of streams, basis is either "Mass" or "Molar", tracer_species is a string, and tracer_concentration is a numpy array
        #generate the total amount of tracer species coming into the system

        total_tracer_in = np.zeros(len(inlet_streams[0].flowrate[0]))
        #a len function in stream would be nice!
        
        for stream in inlet_streams:
            
            
            if basis == "Molar" or basis == "Mass":         #Mass not working yet
                
                total_tracer_in += getattr(stream, "calcSpecies%sFlowrate" % basis)(tracer_species)

                 
        #try to break the error around outlet_tracer_concentration == 0
        outlet_tracer_concentration[outlet_tracer_concentration<=0] = np.nan
  
        total_outlet_flowrate = total_tracer_in/outlet_tracer_concentration
        
        #create the return Stream
        return Stream(name, flowrate = (total_outlet_flowrate,"mol/s"), basis = "molar")

    def generate_inlet_outlet_species_flows(self, name_qualifier = None):
        """Creates inlet and outlet flowrates for the species of interest in the self.proc_species list"""
        if self.inlet_streams ==  [] or self.outlet_streams ==  []:
            raise BadStreamError, "Inlet/outlet streams not yet defined"

        if name_qualifier ==  None:
            outlet_name = 'outlet'
            inlet_name = 'inlet'

        else:
            if type(name_qualifier) !=  str:
                raise lflException, "The new column name qualifier must be a string"
            else:
                outlet_name = 'outlet_%s' % name_qualifier
                inlet_name = 'inlet_%s' % name_qualifier
        
        for specie in self.proc_species:
        
            self['%s_%s' % (specie, inlet_name)] = self.calc_species_molar_feedrates(self.inlet_streams, specie)
            self['%s_%s' % (specie,outlet_name)] = self.calc_species_molar_feedrates(self.outlet_streams, specie)
            self.units['%s_%s' % (specie, inlet_name)] = 'mol/s'
            self.units['%s_%s' % (specie,outlet_name)] = 'mol/s'

    def calc_species_molar_feedrates(self, stream_list, specie):
        """Returns the species molar feedrates in a given set of streams"""
        tot = 0
        for stream in stream_list:
            if not isinstance(stream, Stream):

                raise BadStreamError, "%s is not a stream, which is required for elemental calculations" % stream
            
            tot +=  stream.calcSpeciesMolarFlowrate(specie)
        return tot

    def _calc_conversion(self, inlet_carbon = None, outlet_carbon = None):
        """Calculates a conversion from the given lists of names"""
        if inlet_carbon is None or outlet_carbon is None:
            pass #raise an error here on a problematic list


        #OK, I know, this seems rather unsophisticated and stupid.  However, it is REQUIRED to have simple functions that return an array for uncertainty determination purposes.
        try:
            conv_val = outlet_carbon/inlet_carbon
        except KeyError:
            raise df.NoColumnError, "The desired column is not in the data frame"
        return conv_val

    def _calc_normalized_comp(self,species_to_normalize, excluded_list = None):
        """Calculates a normalized composition given a list of inerts to exclude"""
        try:
            if excluded_list is None:
                return self[species_to_normalize]
            arr_list = []
            for spec in excluded_list:
                arr_list.append(self[spec])
            return self[species_to_normalize]/(1.0-sum(arr_list))
        except KeyError:
            raise df.NoColumnError, "One of the columns in the excluded list or the species to normalize was not defined before attempting normalization"
        except ValueError:
            raise lflException, "You tried to divide by zero!"

    def inlet_enthalpy(self, units):
        """Returns the inlet enthalpy in the process"""
        H = 0.0
        for stream in self.inlet_streams:
            H +=  stream.get_enthalpy(units)
        return (H, units)

    def outlet_enthalpy(self, units):
        """Returns the outlet enthalpy in the process"""
        H = 0.0
        for stream in self.outlet_streams:
            H +=  stream.get_enthalpy(units)
        return (H, units)

    def enthalpy_change(self, units):
        """Returns the enthalpy change in the process"""
        #Create a reactor object with inlet and outlet streams
        r = Reactor(self.inlet_streams, self.outlet_streams)
        return r.enthalpy_change(units)

    def generate_enthalpy_change(self, units):
        """Generates a new column with enthalpy change in the data series"""
        self['delta_H'], self.units['delta_H'] = self.enthalpy_change(self, units)
        

    def collected_enthalpy(self, streams, units):
        """Returns the enthalpy for a collection of streams"""
        H = 0.0
        for stream in streams:
            H +=  stream.get_enthalpy(units)
        return (H, units)

    def enthalpy_diff(self, stream_list_1, stream_list_2, units):
        return self.collected_enthalpy(stream_list_2, units) - self.collected_enthalpy(stream_list_1, units)

class GasifierProcTS(ProcTS):
    """Timeseries data for the gasification chemical process"""
    def __init__(self, start, end, data = None, units_dict = None):
        ProcTS.__init__(self, start, end, data, units_dict)

    

    def generate_carbon_conversions(self, inlet_qualifier = "", qualifier = ""):
        """This calculates the four "standard" conversions/yields that we look for in the biomass gasification on an instantaneous basis"""
        try:
            #CO Yield
            self['CO_yield%s' % qualifier] = self._calc_conversion(self['C_inlet%s' % inlet_qualifier]-self['CO2_inlet%s'% inlet_qualifier],self['CO_outlet%s' % qualifier])
            self.units['CO_yield%s' % qualifier] = None
            #filter these values to throw out the high ones
            self.filter_vals('CO_yield%s' % qualifier, 1.00, 'high')
            self.filter_vals('CO_yield%s' % qualifier, 0.00, 'low')
       
            #Good conversion
            self['X_good%s' % qualifier] = self._calc_conversion(self['C_inlet%s' % inlet_qualifier]-self['CO2_inlet%s' % inlet_qualifier],sum([self['CO_outlet%s' % qualifier],self['CO2_outlet%s' % qualifier]])-self['CO2_inlet%s' % inlet_qualifier])
            self.units['X_good%s' % qualifier] = None
            self.filter_vals('X_good%s' % qualifier, 1.00, 'high')
            self.filter_vals('X_good%s' % qualifier, 0.00, 'low')

            #Standard conversion
            self['X_std%s' % qualifier] = self._calc_conversion(self['C_inlet%s' % inlet_qualifier]-self['CO2_inlet%s' % inlet_qualifier],sum([self['CO_outlet%s' % qualifier],self['CO2_outlet%s' % qualifier],self['CH4_outlet%s' % qualifier]])-self['CO2_inlet%s' % inlet_qualifier])
            self.units['X_std%s' % qualifier] = None
            self.filter_vals('X_std%s' % qualifier, 1.00, 'high')
            self.filter_vals('X_std%s' % qualifier, 0.00, 'low')


            #Total conversion to gas
            self['X_tot%s' % qualifier] = self._calc_conversion(self['C_inlet%s' % inlet_qualifier]-self['CO2_inlet%s' % inlet_qualifier],self['C_outlet%s' % qualifier]-self['CO2_inlet%s' % inlet_qualifier])
            self.units['X_tot%s' % qualifier] = None
            self.filter_vals('X_tot%s' % qualifier, 1.00, 'high')
            self.filter_vals('X_tot%s' % qualifier, 0.00, 'low')

        except KeyError:
            raise ConversionError, "The necessary columns to calculate CO yield have not been generated yet."

    def generate_normalized_compositions(self, name_qualifier = ""):
        """Generate normalized species from the list of process species"""
        norm_flow = sum([self["%s_outlet%s" % (spec, name_qualifier)] for spec in self.proc_species]) - sum([self["%s_outlet%s" % (spec, name_qualifier)] for spec in self.inert_species])

        for species in self.proc_species:
            self['%s_normalized' % species] = self["%s_outlet%s" % (species, "%s" % name_qualifier)]/norm_flow

    def calc_tar_rate(self, exit_stream, name_qualifier = ""):
        """Calculate the level of tar in mg/Nm^3 exiting the gasifier"""
        #I'm just going to start by assuming that the outlet flowrate is molar -- I can make it more general later
	tar_list = ['C6H6', 'C7H8', 'C10H8']
        inclusive_tar_list = ['C2H2', 'C2H4', 'C2H6', 'C3H8', 'C3H6', 'C6H6', 'C7H8', 'C10H8']

        outlet_vol_rate = exit_stream.flowrate[0] * 0.0224	#Nm^3/s, assuming mol/s for basis of original flowrate -- make it more general later

        total_tar = np.zeros(self.nRows)
        total_tar_incl = np.zeros(self.nRows)
        #total tar mass rate
        
        for molecule in tar_list:
            
            total_tar += Molecule(molecule).MW() * self['%s_outlet%s' % (molecule, name_qualifier)]

        for molecule in inclusive_tar_list:
            
            total_tar_incl += Molecule(molecule).MW() * self['%s_outlet%s' % (molecule, name_qualifier)]

        self['tar_loading'] = total_tar/outlet_vol_rate*1000.0
        self['tar_loading_incl'] = total_tar_incl/outlet_vol_rate*1000.0
        self.units['tar_loading'] = 'mg/m^3'
        self.units['tar_loading_incl'] = 'mg/m^3'


    def calc_space_time(self, reactor_vol, excluded_species):
        """Calculates the inlet space time of the reactor based on the inlet streams"""
        conv = uc.UnitConverter()
        vol = conv.convert_units(reactor_vol[0], reactor_vol[1], 'm^3')
        V_dot = 0
        for stream in self.inlet_streams:
            for species in stream.composition:
                if species not in excluded_species:
                    V_dot += stream.calcSpeciesVolumetricFlowrate(species)

        #Create a mixer
        mix = Mixer(inlets = self.inlet_streams)
        #need to get the volumetric flowrate of this stream -- add a function to stream to calculate a gas volume if gas phase -- generalize later

        tau = vol/V_dot
        self['space_time'] = tau
        self.units['space_time'] = 's'

    def calc_min_residence_time(self):
        """Calculates the minimum bound on the residence time, assuming complete conversion and heat up at the instant materials enter the reactor"""
        pass      




class RunInformation:
    """Container class for information about experimental runs -- wrapper class around the dictionary object"""
    def __init__(self):
        self.info = {}
        self.start = None
        self.end = None


    def SQL_load(self, interface, table, run_id):
        """Load the data in from the given table into member objects"""
        query = SQL.select_Query(objects = ['*'], table = table, condition_list = ['run_id=%s' % run_id])
        results = interface.query(query)

        #results will be a list of dicts
        for key, value in results[0].items():
            self.info[key] = value               #we'll use a glossary entry to convert the field names to variable names in the program, so I can choose whatever I want

    def __getitem__(self, index):
        try:
            return self.info[index]
        except KeyError:
            raise lflException, "The desired key is not in the run information container"

    def __setitem__(self, index, value):
        try:
            self.info[index] = value
        except KeyError:
            raise lflException, "The desired key is not in the run information container"
        
if __name__ == '__main__':

    pass
   
