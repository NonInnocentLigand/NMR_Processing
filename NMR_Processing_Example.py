import nmrglue as ng
from matplotlib import pyplot as plt
from matplotlib import colormaps
import numpy as np
import os
import pathlib as pl
import pandas as pd
from NMRInitializer import SingleNMRInitializer, SingleVarian, MultipleNMRInitializer

single_nmrfile_dir = ""  # NMR folder filepath containing fid and procpar files
multiple_nmrfile_dir = ""  # Directory containing only multiple NMR folders

# Make a SingleVarian object to process and plot the data
varian = SingleVarian(single_nmrfile_dir)
varian.to_plot(show_plot=True)

# Make a MultipleNMRInitializer object to process and plot multiple datasets
multiple_varian = MultipleNMRInitializer(multiple_nmrfile_dir)
multiple_varian.to_3d_plot()
# Use the to_plot_concentration_over_time() function to generate a plot of the corrected concentrations for a
# given peak realtive to the internal standard
multiple_varian.to_plot_concentration_over_time(set_peak_width_ppm=0.01, integration_value_int_std_peak=3,
                                                internal_std_peak_ppm=3.78, peak_to_analyze=1.3,
                                                number_of_protons_peak_to_analyze=2, show_plot=True, save_fig=True)
