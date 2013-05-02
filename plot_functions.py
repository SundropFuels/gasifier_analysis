import matplotlib.pyplot as plt

def plot_temperatures(data):
    plt.plot(data['timestamp'], data['temp_steam_reactor_entry'], label='Steam Into Reactor')
    plt.plot(data['timestamp'], data['temp_full_port_body_valve'], label='Valve Body')
    plt.plot(data['timestamp'], data['temp_skin_steam_line_to_lance'], label='Line Above Valve')
    plt.xlabel('Time')
    plt.ylabel('Temperature (C)')
    plt.title('Gasification January 11, 2013')
    plt.legend(loc='best')
    plt.show()
    
def plot_pressures(data):
    plt.plot(data['timestamp'], data['pressure_ash_knockout_vessel'], label='System Pressure')
    plt.plot(data['timestamp'], data['pressure_entrainment'], label='Cross Brush')
    plt.plot(data['timestamp'], data['pressure_feeder_vessel'], label='Down Bed')
    plt.xlabel('Time')
    plt.ylabel('Pressure (psig)')
    plt.title('Gasification January 11, 2013')
    plt.legend(loc='best')
    plt.show()

def plot_feed(data):
    plt.plot(data['timestamp'], data['mass_flow_brush_feeder'], label='Biomass Flow Rate')
    plt.plot(data['timestamp'], data['pressure_feeder_vessel']-data['pressure_ash_knockout_vessel'], label='Pressure Differential')
    plt.plot(data['timestamp'], data['mass_flow_entrainment'], label='Cross Brush Flow')
    plt.xlabel('Time')
    plt.ylabel('Flow Rate (lbs/hr) or dP (psi) or SLPM')
    plt.title('Gasification January 11, 2013')
    plt.legend(loc='best')
    plt.show()
    
def plot_gas_analysis(data):
    plt.plot(data['timestamp'], data['Carbon Monoxide_MS'], label='CO')
    plt.plot(data['timestamp'], data['Hydrogen_MS'], label='H2')
    plt.plot(data['timestamp'], data['Carbon Dioxide_MS'], label='CO2')
    plt.plot(data['timestamp'], data['Methane_MS'], label='CH4')
    plt.xlabel('Time')
    plt.ylabel('% Gas Volume')
    plt.title('Gasification January 11, 2013')
    plt.legend(loc='best')
    plt.show()
