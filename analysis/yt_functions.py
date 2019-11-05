import yt
from yt import YTArray
from yt import YTQuantity
from astropy import constants as const

import numpy as np

#ds = yt.load('../../simulations/isothermal_tctf_1.0/DD0000/DD0000')
ds = yt.load('../../simulations/isocool_tctf_1.0_beta_10.0/DD0000/DD0000')

def _gasentropy(field, data):
    gamma = 5./3.
    mh     = YTQuantity(const.m_p.cgs.value, 'g')
    mu = 1.22
    
    p = data[('gas', 'pressure')]
    n = data[('gas', 'density')] / (mu * mh)
    return  p / np.power(n, gamma)

def _crentropy(field, data):
    crgamma = 4./3.
    mh     = YTQuantity(const.m_p.cgs.value, 'g')
    mu = 1.22
    
    pcr = data[('gas', 'cr_pressure')]
    n = data[('gas', 'density')] / (mu * mh)
    return  pcr / np.power(n, crgamma)

def _log_gas_entropy(field, data):
    return np.log(data[('gas', 'gas_entropy')])

def _log_total_entropy(field, data):
    ge = data[('gas', 'gas_entropy')]
    cre = data[('gas', 'cr_entropy')]
    return np.log(ge) + data[('gas', 'cr_eta')]*np.log(cre)

def _crenergy(field, data):
    cre = ds.arr(data[('enzo', 'CREnergyDensity')] / data[('enzo', 'Density')], 'code_velocity**2')
    return cre

def _crpressure(field, data):
    crgamma = 4./3.
    crp = ds.arr(data[('enzo', 'CREnergyDensity')] * (crgamma - 1.0), 'code_velocity**2 * code_density')
    return crp

def _creta(field, data):
    return data[('gas', 'cr_pressure')] / data[('gas', 'pressure')]
    
def _accel_z(field, data):
    accel_unit = ds.length_unit.value / ds.time_unit.value**2
    accel = data[('enzo', 'External_Acceleration_z')] * accel_unit
    return YTArray(accel, 'cm/s**2')

def _ff_time(field, data):
    z = data[('gas', 'z')]
    return np.sqrt(abs(2.0 * z / data[('gas', 'external_acceleration_z')]))

def _cool_ff_ratio(field, data):
    return data[('gas', 'cooling_time')] / data[('gas', 'free_fall_time')]

def _total_pressure(field, data):
    return data[('gas', 'pressure')] + data[('gas', 'magnetic_pressure')]

def load(output_location):
    ds = yt.load(output_location)
    ds.add_field(('gas', 'external_acceleration_z'), function = _accel_z, \
                display_name = 'External Acceleration Z', units = 'cm/s**2')

    ds.add_field(('gas', 'free_fall_time'), function = _ff_time, \
                display_name = 'Free Fall Time', units = 's')

    ds.add_field(('gas', 'tcool_tff_ratio'), function = _cool_ff_ratio, \
                display_name = 'Cooling Time / Free Fall Time', units = '')
    ds.add_field(('gas', 'total_pressure'), function = _total_pressure, \
                 display_name = 'Total Pressure', units = 'dyne/cm**2')
    ds.add_field(('gas', 'gas_entropy'), function = _gasentropy,\
                 display_name = 'Gas Entropy', units = 'cm**2 * keV')

    ds.add_field(('gas', 'log_total_entropy'), function = _log_gas_entropy, \
                 display_name = 'Log(Total Entropy)', units = '')

    if output_location.__contains__('_cr_'):
        ds.add_field(('gas', 'cr_energy'), function = _crenergy, \
                     display_name = 'Cosmic Ray Energy', units = 'erg/g')
        ds.add_field(('gas', 'cr_pressure'), function = _crpressure, \
                     display_name = 'Cosmic Ray Pressure', units = 'dyne/cm**2')
        ds.add_field(('gas', 'cr_eta'), function = _creta, \
                     display_name = 'P$_c$ / P$_g$', units = '')
        ds.add_field(('gas', 'cr_entropy'), function = _crentropy, \
                     display_name = 'Cosmic Ray Entropy', units = 'cm**2 * dyne')
        ds.add_field(('gas', 'log_total_entropy'), function= _log_total_entropy, \
                     display_name ='Log(Total Entropy)', units = '', force_overwrite = True)
    return ds
