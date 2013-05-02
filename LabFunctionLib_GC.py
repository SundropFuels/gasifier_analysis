"""LabFunctionLib.py
   General data processing and handling functions for experimental analysis
   Chris Perkins
   2012-01-23

   Version history:
   Version 0.1 - First porting over from R libraries

"""

import numpy as np
import db_toolbox as SQL
import dataFrame_v2 as df
import datetime
import unitConversion as uc
import Thermo
import matplotlib.pyplot as plt

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
    def __init__(self, start, end, array_dict = None, units_dict = None):
        df.Dataframe.__init__(self, array_dict = array_dict, units_dict = units_dict)
        if type(start) != datetime.datetime or type(end) != datetime.datetime:
            raise TimeError, "The start and end times for a timeseries must be datetime.datetime objects"
        if start > end:
            raise TimeError, "The start time should happen before the end time"

        self.start = start
        self.end = end
        self.uncertainty_dict = {}  #The uncertainty dictionary -- need to add a function to fill this from an SQL table

    def SQL_load(self, db_interface, table = ""):
        if not isinstance(db_interface, SQL.db_interface):
            raise SQLInterfaceError, "The passed interface to the database is not valid"

        start_condition = "timestamp >= '%s'" % self.start.strftime("%Y-%m-%d %H:%M:%S")
        end_condition = "timestamp <= '%s'" % self.end.strftime("%Y-%m-%d %H:%M:%S")
        try:
            self.SQL_load_data(db_interface, table, conditions = [start_condition, end_condition])
        except SQL.NoDBInUseError:
            raise SQLInterfaceError, "No database selected for the passed interface"
        except SQL.InterfaceNotConnectedError:
            raise SQLInterfaceError, "Passed interface is not connected to a server"
        counter = np.arange(self.numrows())
        self.data['counter'] = counter

        #Need to add the units load to this and to the unittests

    def SQL_uc_load(self, db_interface, table = ""):
        pass

    def interpolate_col(self, x, y):
        """Interpolates between None values on a given column ('y') relative to 'x' and adds the interpolated column as 'Name'_interp"""
        if x not in self.data.keys() or y not in self.data.keys():
            raise InterpArgumentError, "The given column name is not in the data frame -- cannot interpolate"

        interp_col = self.data[y].copy()
        #Scroll forward to find the first non-nan/non-None value
        i = 0
        while (self.data[y][i] == None or self.data[y][i] == np.nan) and i < self.nRows:
            i+=1
        if i == self.nRows:
            #The column was empty -- nothing to do
            return interp_col
        begin = i
        while i < self.nRows:
            i += 1
            if i == self.nRows:
                #Trailing nones/nans
                return interp_col
            elif self.data[y][i] == None or self.data[y][i] == np.nan:
                continue
            elif i == self.nRows:
                #Trailing nones/nans
                return interp_col
            elif i - begin == 1:
                #There were no nans or Nones in the mix, move forward
                begin = i
            else:
                #Now interpolate between the values
                end = i
                for j in range(begin+1,end):
                    interp_col[j] = (self.data[y][end] - self.data[y][begin])/(self.data[x][end] - self.data[x][begin])*(self.data[x][j] - self.data[x][begin]) + self.data[y][begin]
                begin = i
        return interp_col

    def calc_rel_time(self):
        rel_time=[]
        for i in self['timestamp']:
            rel_time.append((i-self['timestamp'][0]).seconds)
        rel_time=np.array(rel_time)
        return rel_time                            

    def rid_nones(self, colname):
        for i in range(len(colname)-1):
            if colname[i-1]==None and colname[i]!=None:
                firstval=i
            if colname[i+1]==None and colname[i]!=None:
                lastval=i
        colname[0:firstval]=colname[firstval]
        colname[lastval+1:]=colname[lastval]
        return(colname)
	
    def interp_ga(self):
        self['rel_time']=self.calc_rel_time()
        for i in self.data:
            if i.endswith('_MS') or i.endswith('_GC'):
                colname=i
                try:
                    interpcol=self.interpolate_col('rel_time', colname)
                    interpcol=self.rid_nones(interpcol)
                    self[colname]=interpcol
                except IndexError: pass 

    def calc_water_ms_corrected(self):
        self.data['water_ms_corrected']=self.data['Water_MS']*(1-self.data['Ar_GC']/100)

    def calc_dry_comp(self, colname):
        self.data[colname+'_dry']=self.data[colname]*(1-self.data['Ar_GC']/100)/(
            1-(self.data['water_ms_corrected']/100))

    def calc_dry_comp_ppm(self, colname):
        self.data[colname+'_dry']=0.0001*self.data[colname]*(1-self.data['Ar_GC']/100)/(
            1-(self.data['water_ms_corrected']/100))

    def calc_all_dry_comp(self, ppmlist=['Benzene_MS', 'Hydrogen Sulfide_MS', 'Napthalene_MS', 'Toluene_MS']):
        gaslist=[]
        for i in self.data:
            if i.endswith('_MS') and i != 'Water_MS' and i !='Argon_MS':
                gaslist.append(i)
        for i in gaslist:
            self.calc_dry_comp(i)
        for i in ppmlist:
            self.calc_dry_comp_ppm(i)
    
    def calc_outlet_dry_flowrate(self):
        """Calculate dry outlet gas flow rate in moles/s"""
        flowrate=self.data['mass_flow_argon_tracer']/(self.data['Ar_GC']/100)/24.13/60
        self.data['outlet_dry_flowrate']=flowrate

    def calc_gas_produced_flowrate(self):
        """Calculate product gas flow rate in moles/s"""
##        gas_produced=self.data['outlet_dry_flowrate']-(self.data['mass_flow_argon_tracer']+self.data[
##            'mass_flow_argon_tracer']+self.data['mass_flow_argon_tracer']*self.data['N2_GC']/self.data['Ar_GC'])/24.134/60
        gas_produced=self.data['outlet_dry_flowrate']-(self.data['mass_flow_entrainment']+self.data[
            'mass_flow_argon_tracer']+self.data['mass_flow_feed_vessel_pressure'])/24.13/60
        self.data['gas_produced_flowrate']=gas_produced

    def calc_product_flowrates(self):
        """Calculate product gas flow rates in moles/s"""
        gaslist=[]
        for i in self.data:
            if i.endswith('_MS_dry'):
                gaslist.append(i)
        for i in gaslist:
            product_flowrate=self.data[i]/100*self.data['outlet_dry_flowrate']
            self.data[i+'_flowrate']=product_flowrate

    def calc_carbon_in(self, biomass_C=0.52):
        """Calculate inlet carbon flow rate in biomass in moles/s"""
        carbon_in=(self.data['mass_flow_biomass_from_ktron']+0.03)*biomass_C*0.125997881/12.011
        self.data['carbon_in']=carbon_in

    def calc_carbon_out(self):
        c_content={'Methane':1, 'Acetylene':2, 'Carbon Monoxide':1, 'Ethylene':2, 'Ethane':2, 'Propylene':3,
                   'Propane':3, 'Carbon Dioxide':1, 'Benzene':6, 'Toluene':7, 'Napthalene':10}
        for k in c_content:
            carbon_out=self.data[k+'_MS_dry_flowrate']*c_content[k]
            self.data[k+'_carbon_out']=carbon_out
        carbon_out_total=range(len(self.data['timestamp']))
        for i in carbon_out_total:
            carbon_out_total[i]=0
            carbon_out_total=np.array(carbon_out_total)
        for i in self.data:
            if i.endswith('_carbon_out'):
                carbon_out_total=carbon_out_total+self.data[i]
        self.data['carbon_out_total']=carbon_out_total

    def calc_volumetric_flow_rate(self):
        """Calculates actual flow rate in L/s"""
        pass

    def calc_dead_time(self):
        """Calculates dead time in seconds from experimental RTD tests"""
        pass

    def calc_average_feed_rate(self):
        return np.average(self.data['mass_flow_biomass_from_ktron'])

    def calc_std_feed_rate(self):
        return np.std(self.data['mass_flow_biomass_from_ktron'])
        
    def calc_total_c_in(self):
        return np.sum(self.data['carbon_in'])

    def calc_total_c_out(self):
        return np.sum(self.data['carbon_out_total'])

    def calc_carbon_x(self):
        return np.sum(self.data['carbon_out_total'])/np.sum(self.data['carbon_in'])

    def calc_y_co(self):
        return np.sum(self.data['Carbon Monoxide_MS_dry_flowrate'])/np.sum(self.data['carbon_in'])

    def calc_ratio_co_h2(self):
        return np.sum(self.data['Carbon Monoxide_MS_dry_flowrate'])/np.sum(self.data['Hydrogen_MS_dry_flowrate'])

    def calc_ratio_co_co2(self):
        return np.sum(self.data['Carbon Monoxide_MS_dry_flowrate'])/np.sum(self.data['Carbon Dioxide_MS_dry_flowrate'])

    def calc_ratio_co_ch4(self):
        return np.sum(self.data['Carbon Monoxide_MS_dry_flowrate'])/np.sum(self.data['Methane_MS_dry_flowrate'])
        
    def convert_col_units(self, column, units):
        """Will convert the given column into the new units from the existing units"""
        
        if column not in self.data.keys():
            raise df.NoColumnError, "The requested column is not in the data frame."

        if self.units == {}:
            raise NoUnitsError, "The units have not been defined for this column"

        conv = uc.UnitConverter()

        try:
            self.data[column] = conv.convert_units(self.data[column], self.units[column], units)
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
        
        
class Stream:
    species = {}
    species['CO']={'C':1,'O':1}
    species['CO2'] = {'C':1,'O':2}
    species['CH4'] = {'C':1, 'H':4}
    species['H2'] = {'H':2}
    species['C2H2'] = {'C':2,'H':2}
    species['C2H4'] = {'C':2,'H':4}
    species['C2H6'] = {'C':2,'H':6}
    species['H2O'] = {'H':2,'O':1}
    species['AR'] = {'AR':1}
    species['N2'] = {'N':2}
    


    def __init__(self, name, flowrate = None, composition = None, basis = "molar", temperature = None, pressure = None, density=None, compressible = None):

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

        self.special_species = {}
        

    def elementalFactor(self, element):
        """Returns the units of element/basis unit of a given element in the feed"""
        factor = 1.0
        for (specie, fraction) in self.composition.items():
            try:
                spe_dict = Stream.species[specie]
            except KeyError:
                try:
                    spe_dict = self.special_species[specie]
                except KeyError:
                    raise SpeciesNotDefinedError, "%s does not have an elemental breakdown definition...check special_species{}" % specie
            try:
                factor += fraction*spe_dict[element]
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
        if self.basis == "molar":
            #First need to convert to mol/s
            conv = uc.UnitConverter()
            val = conv.convert_units(self.flowrate[0], self.flowrate[1], 'mol/s')
            return val*self.elementalFactor(element)
        elif self.basis == "mass":
            #First need to convert to g/s
            conv = uc.UnitConverter()
            val = conv.convert_units(self.flowrate[0], self.flowrate[1], 'g/s')
            return val*self.elementalFactor(element)/Element.MW[element]
        elif self.basis == "gas_volume":
            #First need to convert to m^3/s
            conv = uc.UnitConverter()
            val = conv.convert_units(self.flowrate[0], self.flowrate[1], 'm^3/s')
            p = conv.convert_units(self.pressure[0], self.pressure[1], 'kg/s^2/m')
            T = conv.convert_units(self.temperature[0], self.temperature[1], 'K')
            #Need to convert the temperature to K
            return val*p/(8.314*T)*self.elementalFactor(element)
        elif self.basis == "liquid_volume":
            #First need to convert the density to the same units as the liquid
            pass
        else:
            pass

    def calcSpeciesMolarFlowrate(self, species):
        """Calculates the molar feed rate (mol/s) of a species for a given stream flowrate"""
        if species not in self.composition.keys():
            return np.zeros(len(self.flowrate[0]))

        if self.basis == "molar":
            conv = uc.UnitConverter()
            val = conv.convert_units(self.flowrate[0], self.flowrate[1], 'mol/s')
            return val*self.composition[species]
        elif self.basis == "mass":
            #Need to convert to molar amounts based on the molecular weight
            conv = uc.UnitConverter()
            val = conv.convert_units(self.flowrate[0], self.flowrate[1], 'g/s')
            
            if species in SpecialMolecule.MW.keys():
                MW = SpecialMolecule.MW[species]
            else:
                MW = 0
                breakdown = _parse_species(species)
                try:
                    for ele, v in breakdown.items():
                        MW += v*Element.MW[ele]
                except KeyError:
                    raise BadStreamError, "%s does not have an entry in the Element molecular weight dictionary" % ele
            return val*self.composition[species]/MW
        elif self.basis == "gas_volume":
            conv = uc.UnitConverter()
            val = conv.convert_units(self.flowrate[0], self.flowrate[1], 'm^3/s')
            p = conv.convert_units(self.pressure[0], self.pressure[1], 'kg/s^2/m')
            T = conv.convert_units(self.temperature[0], self.temperature[1], 'K')
            return val*self.composition[species]*p/(8.314*T)
    
    
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
                if i+1 == len(species_str) or species_str[i+1].isupper():
                    #Allows for single character names, like CH4
                    current_element = "".join([current_element, current_char])
                    if current_element in parsed_elements.keys():
                        parsed_elements[current_element] += 1
                    else:
                        parsed_elements[current_element] = 1
                    current_element = ""
                
                
                else:
                
                    current_element = "".join([current_element, current_char])
                i += 1
                continue

            elif current_char.isdigit():
                #we have gotten to the end of an element name
                print current_element
                while i < len(species_str) and species_str[i].isdigit():
                    print current_char
                    current_char = species_str[i]
                    current_number_str = "".join([current_number_str, current_char])
                    i += 1
                if current_number_str == '':
                    raise BadCharacterError, "Each element must have a number associated with it"
                if current_element in parsed_elements.keys():
                    parsed_elements[current_element] += int(current_number_str)
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
            raise UnitConversionError, "The provided units must be a strong for the entropy function"
        self._calc_entropy()
        conv = uc.UnitConverter()
        return conv.convert_units(self.entropy[0], self.entropy[1], units)


    def _calc_enthalpy(self):
        """Calculates the stream enthalpy and stores it in self.enthalpy"""
        tg = Thermo.ThermoGenerator()
       
        H = 0.0
        #Calculate the enthalpy of each stream component and add them all together

        for sp in self.composition.keys():
            if self.basis == "molar":
                conv = uc.UnitConverter()
                val = conv.convert_units(self.flowrate[0], self.flowrate[1], 'mol/s')
                val *= self.composition[sp]
                H += val*tg.calc_enthalpy(sp, self.temperature, 'J/mol')
                



            elif self.basis == "mass":
                conv = uc.UnitConverter()
                val = conv.convert_units(self.flowrate[0], self.flowrate[1], 'g/s')
                if species in SpecialMolecule.MW.keys():
                    MW = SpecialMolecule.MW[sp]
                else:
                    MW = 0
                    breakdown = _parse_species(species)
                    try:
                        for ele, v in breakdown.items():
                            MW += v*Element.MW[ele]
                    except KeyError:
                        raise BadStreamError, "%s does not have an entry in the Element molecular weight dictionary" % ele
                val *= self.composition[sp]/MW #This is now the flowrate of this species in mol/s
                H += val*tg.calc_enthalpy(sp, self.temperature, 'J/mol')
               


            elif self.basis == "gas_volume":
                
                conv = uc.UnitConverter()
                val = conv.convert_units(self.flowrate[0], self.flowrate[1], 'm^3/s')
                p = conv.convert_units(self.pressure[0], self.pressure[1], 'kg/s^2/m')
                T = conv.convert_units(self.temperature[0], self.temperature[1], 'K')
                val *= self.composition[sp]*p/(8.314*T)
                H += val*tg.calc_enthalpy(sp, self.temperature, 'J/mol')
                


        self.enthalpy = (H, 'J/s')
      

    def _calc_entropy(self):
        """Calculates the stream entropy and stores it in self.entropy"""
        tg = Thermo.ThermoGenerator()
       
        S = 0.0
        #Calculate the enthalpy of each stream component and add them all together

        for sp in self.composition.keys():
            if self.basis == "molar":
                conv = uc.UnitConverter()
                val = conv.convert_units(self.flowrate[0], self.flowrate[1], 'mol/s')
                val *= self.composition[sp]
                S += val*tg.calc_entropy(sp, self.temperature, 'J/mol')
                #Need to adjust if compressible is true for the stream -- FIX FOR ENTROPY
                if self.compressible == True:
                    p = conv.convert_units(self.pressure[0], self.pressure[1], 'Pa')
                    T = conv.convert_units(self.temperature[0], self.temperature[1], 'K')
                    S += 8.314 * val * np.log(p/101325.0)



            elif self.basis == "mass":
                conv = uc.UnitConverter()
                val = conv.convert_units(self.flowrate[0], self.flowrate[1], 'g/s')
                if species in SpecialMolecule.MW.keys():
                    MW = SpecialMolecule.MW[sp]
                else:
                    MW = 0
                    breakdown = _parse_species(species)
                    try:
                        for ele, v in breakdown.items():
                            MW += v*Element.MW[ele]
                    except KeyError:
                        raise BadStreamError, "%s does not have an entry in the Element molecular weight dictionary" % ele
                val *= self.composition[sp]/MW #This is now the flowrate of this species in mol/s
                S += val*tg.calc_entropy(sp, self.temperature, 'J/mol')
                if self.compressible == True:
                    p = conv.convert_units(self.pressure[0], self.pressure[1], 'Pa')
                    T = conv.convert_units(self.temperature[0], self.temperature[1], 'K')
                    S += 8.314 * val * np.log(p/101325.0)


            elif self.basis == "gas_volume":
                
                conv = uc.UnitConverter()
                val = conv.convert_units(self.flowrate[0], self.flowrate[1], 'm^3/s')
                p = conv.convert_units(self.pressure[0], self.pressure[1], 'kg/s^2/m')
                T = conv.convert_units(self.temperature[0], self.temperature[1], 'K')
                val *= self.composition[sp]*p/(8.314*T)
                S += val*tg.calc_entropy(sp, self.temperature, 'J/mol')
                S += 8.314 * val * np.log(p/101325.0)
        
        self.entropy = (S, 'J/mol/K')


class PhysicalConstants:
    """Class to hold global physical constant variables"""
    std_temp = 298   #In K
    std_pressure = 100000 #In Pa

class Element:
    """A physical element.  Will have properties of the elements, plus a handy reference library"""
    Elements = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'AR', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au','Hg', 'Tl', 'Pb','Bi', 'Po', 'At', 'Rn', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu']
    MW = {}
    MW['C'] = 12.00
    MW['H'] = 1.00
    MW['O'] = 16.00
    MW['N'] = 14.00

class Molecule:
    pass

class SpecialMolecule(Molecule):
    """This class holds specially designated molcules, so that I can calculate molecular weights without too much of a problem"""
    MW = {}
    MW['LigH'] = 0
    MW['LigO'] = 0
    MW['LigC'] = 0
    MW['Cell'] = 162
    MW['HemiC'] = 0

class ProcTS(ts_data):

    """Timeseries data for chemical processes"""
    def __init__(self, start, end, array_dict = None, units_dict = None):
        ts_data.__init__(self, start, end, array_dict, units_dict)
        
        self.inlet_streams = []
        self.outlet_streams = []
        self.proc_elements = []
        self.proc_species = []
        self.inert_species = []

    def generate_inlet_outlet_elemental_flows(self, name_qualifier = None):
        """Generates the standard inlet and outlet elemental flows for elements in self.proc_elements"""
        if self.inlet_streams == [] or self.outlet_streams == []:
            raise BadStreamError, "Inlet/outlet streams not yet defined"
        
        if name_qualifier == None:
            outlet_name = 'outlet'
            inlet_name = 'inlet'

        else:
            if type(name_qualifier) != str:
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
                
            tot += stream.calcElementalMolarFlowrate(element)

        return tot

    def generate_inlet_outlet_species_flows(self, name_qualifier = None):
        """Creates inlet and outlet flowrates for the species of interest in the self.proc_species list"""
        if self.inlet_streams == [] or self.outlet_streams == []:
            raise BadStreamError, "Inlet/outlet streams not yet defined"

        if name_qualifier == None:
            outlet_name = 'outlet'
            inlet_name = 'inlet'

        else:
            if type(name_qualifier) != str:
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

        tot += stream.calcSpeciesMolarFlowrate(specie)
        return tot

    def _calc_conversion(self, inlet_carbon = None, outlet_carbon = None):
        """Calculates a conversion from the given lists of names"""
        if inlet_carbon == None or outlet_carbon == None:
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
            return self[species_to_normalize]/sum(arr_list)
        except KeyError:
            raise df.NoColumnError, "One of the columns in the excluded list or the species to normalize was not defined before attempting normalization"
        except ValueError:
            raise lflException, "You tried to divide by zero!"

    def inlet_enthalpy(self, units):
        """Returns the inlet enthalpy in the process"""
        H = 0.0
        for stream in self.inlet_streams:
            H += stream.get_enthalpy(units)
        return (H, units)

    def outlet_enthalpy(self, units):
        """Returns the outlet enthalpy in the process"""
        H = 0.0
        for stream in self.outlet_streams:
            H += stream.get_enthalpy(units)
        return (H, units)

    def enthalpy_change(self, units):
        """Returns the enthalpy change in the process"""
        H = self.outlet_enthalpy(units) - self.inlet_enthalpy(units)
        return (H, units)

    def generate_enthalpy_change(self, units):
        """Generates a new column with enthalpy change in the data sereies"""
        self['delta_H'], self.units['delta_H'] = self.enthalpy_change(self, units)
        

    def collected_enthalpy(self, streams, units):
        """Returns the enthalpy for a collection of streams"""
        H = 0.0
        for stream in streams:
            H += stream.get_enthalpy(units)
        return (H, units)

    def enthalpy_diff(self, stream_list_1, stream_list_2, units):
        return self.collected_enthalpy(stream_list_2, units) - self.collected_enthalpy(stream_list_1, units)

class GasifierProcTS(ProcTS):
    """Timeseries data for the gasification chemical process"""
    def __init__(self, start, end, array_dict = {}, units_dict = {}):
        ProcTS.__init__(self, start, end, array_dict, units_dict)

    def generate_carbon_conversions(self):
        """This calculates the four "standard" conversions/yields that we look for in the biomass gasification on an instantaneous basis"""
        try:
            #CO Yield
            self['CO_yield'] = self._calc_conversion(self['C_inlet'],self['CO_outlet'])
            self.units['CO_yield'] = None
        
       
            #Good conversion
            self['X_good'] = self._calc_conversion(self['C_inlet'],sum([self['CO_outlet'],self['CO2_outlet']]))
            self.units['X_good'] = None

            #Standard conversion
            self['X_std'] = self._calc_conversion(self['C_inlet'],sum([self['CO_outlet'],self['CO2_outlet'],self['CH4_outlet']]))
            self.units['X_std'] = None

            #Total conversion to gas
            self['X_tot'] = self._calc_conversion(self['C_inlet'],self['C_outlet'])
            self.units['X_tot'] = None

        except df.NoColumnError:
            raise ConversionError, "The necessary columns to calculate CO yield have not been generated yet."

    def generate_normalized_compositions(self):
        """Generate normalized species from the list of process species"""
        for species in self.proc_species:
            self['%s_normalized' % species] = self._calc_normalized_comp(species_to_normalize = species, excluded_list = self.inert_species)


class RunInformation:
    """Container class for information about experimental runs -- wrapper class around the dictionary object"""
    def __init__(self):
        self.info = {}
        self.start = None
        self.end = None


    def SQL_load(self, interface, table, run_id):
        """Load the data in from the given table into member objects"""
        query = SQL.select_Query(objects = ['*'], table = table, condition_list = [run_id])
        results = interface.query(query)

        #results will be a list of dicts
        for key, value in results[0].items():
            self.info[key] = value               #we'll use a glossary entry to convert the field names to variable names in the program, so I can choose whatever I want

        
if __name__=='__main__':

    import collections
    table=collections.OrderedDict({})
    
    d=collections.OrderedDict({})
##    d['GAS20121113a']=['2012-11-20 3:25:00', '2012-11-20 3:29:30']
##    d['GAS20121113e']=['2012-11-20 5:40:00', '2012-11-20 5:46:00']
##    d['GAS20121113f']=['2012-11-27 11:53:00', '2012-11-27 12:03:00']
##    d['GAS20121113g']=['2012-11-21 16:53:00', '2012-11-21 17:03:00']
##    d['GAS20121113c']=['2012-11-16 23:31:00', '2012-11-16 23:36:00']
##    d['GAS20121113d']=['2012-11-16 20:42:00', '2012-11-16 20:49:00']
##    d['GAS20121113h']=['2012-11-20 9:48:00', '2012-11-20 9:53:00']
##    d['GAS20121113i']=['2012-11-21 13:32:00', '2012-11-21 13:37:00']
##    d['GAS20121113j']=['2012-11-21 22:11:00', '2012-11-21 22:19:00']
##    d['GAS20121113k']=['2012-11-27 16:15:00', '2012-11-27 16:29:00']
##    d['GAS20121113c']=['2012-12-7 13:44:00', '2012-12-7 14:46:00']

    d['GAS20121113c']=['2012-12-17 16:45:00', '2012-12-17 16:50:00']    
    for k in d:
        start=d[k][0]
        end=d[k][1]
        dbconn=SQL.db_interface(host='io', user='jmp_user', passwd='jmpyourface')
        dbconn.connect()
        dtstart=datetime.datetime.strptime(start, '%Y-%m-%d %H:%M:%S')
        dtend=datetime.datetime.strptime(end, '%Y-%m-%d %H:%M:%S')
        print k
        print 'Start: ' + d[k][0] + '\n' + 'End: ' + d[k][1]
        data=ts_data(dtstart, dtend)
        data.SQL_load(dbconn, table='gasifier_data_HR_view')
##        data.interp_ga()
##        data.calc_water_ms_corrected()
##        data.calc_all_dry_comp()
##        data.calc_outlet_dry_flowrate()
##        data.calc_gas_produced_flowrate()
##        data.calc_product_flowrates()
##        data.calc_carbon_in()
##        data.calc_carbon_out()
##        table[k]=collections.OrderedDict({})
##        table[k]['Average Feed Rate']=round(data.calc_average_feed_rate(),2)
##        table[k]['StD Feed Rate']=round(data.calc_std_feed_rate(),2)
##        table[k]['Carbon In']=round(data.calc_total_c_in(),1)
##        table[k]['Carbon Out']=round(data.calc_total_c_out(),1)
##        table[k]['Carbon Conversion']=round(data.calc_carbon_x(),3)
##        table[k]['Carbon Yield']=round(data.calc_y_co(),3)
##        table[k]['CO/H2']=round(data.calc_ratio_co_h2(),3)
##        table[k]['CO/CO2']=round(data.calc_ratio_co_co2(),3)
##        table[k]['CO/CH4']=round(data.calc_ratio_co_ch4(),3)
##        print 'Average Feed Rate: ' + str(round(np.average(data['mass_flow_biomass_from_ktron']),2)) + ' lbs/hr'
##        print 'StD Feed Rate: ' + str(round(np.std(data['mass_flow_biomass_from_ktron']),2)) + ' lbs/hr'
##        print 'Carbon in: ' +str(round(np.sum(data['carbon_in']),1)) + ' moles'
##        print 'Carbon out: ' +str(round(np.sum(data['carbon_out_total']),1)) + ' moles'
##        print 'Total C X: ' +str(round(np.sum(data['carbon_out_total'])/np.sum(data['carbon_in']),3))
##        print 'Y CO: ' + str(round(np.sum(data['Carbon Monoxide_MS_dry_flowrate'])/np.sum(data['carbon_in']),3))
##        print 'CO/H2: ' + str(round(np.sum(data['Carbon Monoxide_MS_dry_flowrate'])/np.sum(data['Hydrogen_MS_dry_flowrate']),3))
##        print 'CO/CO2: ' + str(round(np.sum(data['Carbon Monoxide_MS_dry_flowrate'])/np.sum(data['Carbon Dioxide_MS_dry_flowrate']),3))
##        print 'CO/CH4: ' + str(round(np.sum(data['Carbon Monoxide_MS_dry_flowrate'])/np.sum(data['Methane_MS_dry_flowrate']),3)) +'\n'
##        plt.plot(data['timestamp'], data['mass_flow_biomass_from_ktron'])
##        plt.plot(data['timestamp'], data['setpoint_biomass_feedrate'])
##        plt.plot(data['timestamp'], data['Carbon Monoxide_MS']/10)
##        plt.plot(data['timestamp'], data['carbon_out_total'])
##        plt.show()

##    with open('analysis.csv', 'w') as f:
##        f.write('Run ID,')
##        for i in table[k]:
##            f.write(i+',')
##        f.write('\n')
##        for i in table:
##            f.write(i+',')
##            for j in table[i]:
##                f.write(str(table[i][j])+',')
##            f.write('\n')
##    
   
