import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
from matplotlib import animation
from matplotlib.animation import PillowWriter

from scipy.interpolate import interp1d


def constant_time_transform(data, times_in, times_new):
    data_out = []

    for data_set in data:
        data_interp = interp1d(times_in, data_set)
        data_constant_time_points = data_interp(times_new)
        data_out.append(data_constant_time_points)

    return data_out

