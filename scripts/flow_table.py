import pandas as pd
from scipy import interpolate
import numpy as np


def print_flow_table(given_variable, value, flow_type='both'):
    gas_table = pd.read_csv('gas_table.csv')

    mach_sub = np.array(gas_table['Ma'][:100]).astype('float')
    mach_sup = np.array(gas_table['Ma'][101:]).astype('float')

    pres_sub = np.array(gas_table['P/P_0'][:100]).astype('float')
    pres_sup = np.array(gas_table['P/P_0'][101:]).astype('float')

    temp_sub = np.array(gas_table['T/T_0'][:100]).astype('float')
    temp_sup = np.array(gas_table['T/T_0'][101:]).astype('float')

    area_sub = np.array(gas_table['A/A*'][:100]).astype('float')
    area_sup = np.array(gas_table['A/A*'][101:]).astype('float')

    if given_variable == 'M':
        mach = value

        if mach > 1:
            pres = interpolate.interp1d(mach_sub, pres_sub)(value)
            temp = interpolate.interp1d(mach_sub, temp_sub)(value)
            area = interpolate.interp1d(mach_sub, area_sub)(value)
        else:
            pres = interpolate.interp1d(mach_sup, pres_sup)(value)
            temp = interpolate.interp1d(mach_sup, temp_sup)(value)
            area = interpolate.interp1d(mach_sup, area_sup)(value)

    elif given_variable == 'P':
        pres = value

        if flow_type == 'subsonic':
            mach = interpolate.interp1d(pres_sub, mach_sub)(value)
            temp = interpolate.interp1d(pres_sub, temp_sub)(value)
            area = interpolate.interp1d(pres_sub, area_sub)(value)
        elif flow_type == 'supersonic':
            mach = interpolate.interp1d(pres_sup, mach_sup)(value)
            temp = interpolate.interp1d(pres_sup, temp_sup)(value)
            area = interpolate.interp1d(pres_sup, area_sup)(value)
        else:
            pres = [value, value]
            mach = [interpolate.interp1d(pres_sub, mach_sub)(value), interpolate.interp1d(pres_sup, mach_sup)(value)]
            temp = [interpolate.interp1d(pres_sub, temp_sub)(value), interpolate.interp1d(pres_sup, temp_sup)(value)]
            area = [interpolate.interp1d(pres_sub, area_sub)(value), interpolate.interp1d(pres_sup, area_sup)(value)]

    elif given_variable == 'T':
        temp = value

        if flow_type == 'subsonic':
            mach = interpolate.interp1d(temp_sub, mach_sub)(value)
            pres = interpolate.interp1d(temp_sub, pres_sub)(value)
            area = interpolate.interp1d(temp_sub, area_sub)(value)
        elif flow_type == 'supersonic':
            mach = interpolate.interp1d(temp_sup, mach_sup)(value)
            pres = interpolate.interp1d(temp_sup, pres_sup)(value)
            area = interpolate.interp1d(temp_sup, area_sup)(value)
        else:
            temp = [value, value]
            mach = [interpolate.interp1d(temp_sub, mach_sub)(value), interpolate.interp1d(temp_sup, mach_sup)(value)]
            pres = [interpolate.interp1d(temp_sub, pres_sub)(value), interpolate.interp1d(temp_sup, pres_sup)(value)]
            area = [interpolate.interp1d(temp_sub, area_sub)(value), interpolate.interp1d(temp_sup, area_sup)(value)]

    elif given_variable == 'A':
        area = value

        if flow_type == 'subsonic':
            mach = interpolate.interp1d(area_sub, mach_sub)(value)
            pres = interpolate.interp1d(area_sub, pres_sub)(value)
            temp = interpolate.interp1d(area_sub, temp_sub)(value)
        elif flow_type == 'supersonic':
            mach = interpolate.interp1d(area_sup, mach_sup)(value)
            pres = interpolate.interp1d(area_sup, pres_sup)(value)
            temp = interpolate.interp1d(area_sup, temp_sup)(value)
        else:
            area = [value, value]
            mach = [interpolate.interp1d(area_sub, mach_sub)(value), interpolate.interp1d(area_sup, mach_sup)(value)]
            pres = [interpolate.interp1d(area_sub, pres_sub)(value), interpolate.interp1d(area_sup, pres_sup)(value)]
            temp = [interpolate.interp1d(area_sub, temp_sub)(value), interpolate.interp1d(area_sup, temp_sup)(value)]

    if flow_type != 'both' or given_variable == 'M':
        print(f'Mach: {mach} \nPressure Ratio: {pres} \nTemperature Ratio: {temp} \nArea Ratio: {area}')
    else:
        print(f'Subsonic Solution\nMach: {mach[0]} \nPressure Ratio: {pres[0]} \nTemperature Ratio: {temp[0]} \nArea Ratio: {area[0]}')
        print(f'\nSupersonic Solution\nMach: {mach[1]} \nPressure Ratio: {pres[1]} \nTemperature Ratio: {temp[1]} \nArea Ratio: {area[1]}')


def get_flow_values(given_variable, value, flow_type='supersonic'):
    gas_table = pd.read_csv('gas_table.csv')

    mach_sub = np.array(gas_table['Ma'][:100]).astype('float')
    mach_sup = np.array(gas_table['Ma'][101:]).astype('float')

    pres_sub = np.array(gas_table['P/P_0'][:100]).astype('float')
    pres_sup = np.array(gas_table['P/P_0'][101:]).astype('float')

    temp_sub = np.array(gas_table['T/T_0'][:100]).astype('float')
    temp_sup = np.array(gas_table['T/T_0'][101:]).astype('float')

    area_sub = np.array(gas_table['A/A*'][:100]).astype('float')
    area_sup = np.array(gas_table['A/A*'][101:]).astype('float')

    if given_variable == 'M':
        mach = np.array(value)

        if mach > 1:
            pres = interpolate.interp1d(mach_sub, pres_sub)(value)
            temp = interpolate.interp1d(mach_sub, temp_sub)(value)
            area = interpolate.interp1d(mach_sub, area_sub)(value)
        else:
            pres = interpolate.interp1d(mach_sup, pres_sup)(value)
            temp = interpolate.interp1d(mach_sup, temp_sup)(value)
            area = interpolate.interp1d(mach_sup, area_sup)(value)

    elif given_variable == 'P':
        pres = np.array(value)

        if flow_type == 'subsonic':
            mach = interpolate.interp1d(pres_sub, mach_sub)(value)
            temp = interpolate.interp1d(pres_sub, temp_sub)(value)
            area = interpolate.interp1d(pres_sub, area_sub)(value)
        elif flow_type == 'supersonic':
            mach = interpolate.interp1d(pres_sup, mach_sup)(value)
            temp = interpolate.interp1d(pres_sup, temp_sup)(value)
            area = interpolate.interp1d(pres_sup, area_sup)(value)

    elif given_variable == 'T':
        temp = np.array(value)

        if flow_type == 'subsonic':
            mach = interpolate.interp1d(temp_sub, mach_sub)(value)
            pres = interpolate.interp1d(temp_sub, pres_sub)(value)
            area = interpolate.interp1d(temp_sub, area_sub)(value)
        elif flow_type == 'supersonic':
            mach = interpolate.interp1d(temp_sup, mach_sup)(value)
            pres = interpolate.interp1d(temp_sup, pres_sup)(value)
            area = interpolate.interp1d(temp_sup, area_sup)(value)

    elif given_variable == 'A':
        area = np.array(value)

        if flow_type == 'subsonic':
            mach = interpolate.interp1d(area_sub, mach_sub)(value)
            pres = interpolate.interp1d(area_sub, pres_sub)(value)
            temp = interpolate.interp1d(area_sub, temp_sub)(value)
        elif flow_type == 'supersonic':
            mach = interpolate.interp1d(area_sup, mach_sup)(value)
            pres = interpolate.interp1d(area_sup, pres_sup)(value)
            temp = interpolate.interp1d(area_sup, temp_sup)(value)

    return [mach.item(), pres.item(), temp.item(), area.item()]