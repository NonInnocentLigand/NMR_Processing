import nmrglue as ng
from matplotlib import pyplot as plt
from matplotlib import colormaps
import numpy as np
import os
import pathlib as pl
import pandas as pd


class SingleNMRInitializer:
    def __init__(self, file_path):
        self.folder = file_path
        self.dic = 0
        self.data = 0
        self.processed_data = 0
        self.data_ppm = 0
        self.uc = 0  # the unit conversion object
        self.signal_peaks = []  # list of tuples
        # (Peak locations as index of data, Cluster numbers for peaks, Estimated peak scales/linewidths, Estimated peak amplitudes)
        self.integrations = {}  #
        self.normalized_integrations = {}
        self.to_integrations()  # could manually call this method if you want to save computation on initializing instances

    def to_pipe_conversion(self):
        raise NotImplementedError("to_conversion")

    def to_signal_processing(self):
        """take converted data and process it"""
        dic, data = self.to_pipe_conversion()
        self.dic = dic
        self.data = data

        def process(dic, data):
            dic, data = ng.pipe_proc.sp(dic, data, off=0.35, end=0.98, pow=2, c=1.0)
            dic, data = ng.pipe_proc.zf(dic, data, auto=True)
            dic, data = ng.pipe_proc.ft(dic, data, auto=True)
            data = ng.proc_autophase.autops(data, "acme")
            dic, data = ng.pipe_proc.di(dic, data)
            return dic, data

        dic, processed_data = process(dic, data)
        self.processed_data = processed_data
        return processed_data

    def to_ppm(self):
        if not self.dic or self.data.size == 0:  # assumes data is a np array. This could be buggy
            raise NotImplementedError("run to_signal_processing method to get get dic and data")
        uc = ng.pipe.make_uc(self.dic, self.processed_data)
        self.uc = uc
        data_ppm = uc.ppm_scale()

        # This function handles a weird issue with the ranges of ppm produced
        def adjust_ppm_scale(unadjusted_data_ppm: np.array):
            adjusted_data_ppm = unadjusted_data_ppm - unadjusted_data_ppm[unadjusted_data_ppm.shape[0]-1]
            adjusted_data_ppm *= 16
            adjusted_data_ppm -= 2
            return adjusted_data_ppm

        data_ppm = adjust_ppm_scale(data_ppm)
        self.data_ppm = data_ppm
        return data_ppm

    def initialize(self):
        self.to_signal_processing()
        self.to_ppm()

    def to_plot(self):
        raise NotImplementedError("to_plot")

    def to_signal_peaks(self):
        self.initialize()
        signal_peaks = ng.peakpick.pick(data=self.processed_data, pthres=1000000, msep=0.001)
        self.signal_peaks = signal_peaks
        return signal_peaks

    def to_integrations(self):
        self.to_signal_peaks()
        integrations = {}
        for i in self.signal_peaks:
            picked_peak_index = int(i[0])  # gives index of peak in self.data
            lower_bound = picked_peak_index - 20
            upper_bound = picked_peak_index + 20
            peak = self.processed_data[lower_bound:upper_bound + 1]
            integrations[picked_peak_index] = peak
        self.integrations = integrations
        integration_sums = []
        for key in integrations:
            integration_sums.append(integrations[key].sum())
        normalized_integration_sums = integration_sums / min(integration_sums)
        normalized_integrations = dict(zip(list(integrations.keys()), normalized_integration_sums.round(decimals=3)))
        self.normalized_integrations = normalized_integrations

        return integrations, normalized_integrations


class SingleVarian(SingleNMRInitializer):
    def __init__(self, file_path):
        super().__init__(file_path)

    def to_pipe_conversion(self):
        vdic, vdata = ng.varian.read(self.folder, f'{self.folder}/fid', f'{self.folder}/procpar')
        C = ng.convert.converter()
        C.from_varian(vdic, vdata)
        pdic, pdata = C.to_pipe()
        return pdic, pdata

    def to_plot(self, show_plot=False):
        self.initialize()
        plt.plot(self.data_ppm, self.processed_data)
        plt.ylabel("intensity")
        plt.xlabel("σ ppm")
        title = self.folder.split("/")[-1]
        plt.suptitle(f"{title}")
        if show_plot:
            plt.show()

    def to_plot_integrations(self, show_plot=False):
        self.initialize()
        integrations, normalized_integrations = self.to_integrations()
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(self.data_ppm, self.processed_data)

        for k in integrations:
            peak = integrations[k]
            peak_sum = integrations[k].sum()
            lower_bound = k - 20
            upper_bound = k + 20
            x_axis = self.data_ppm[lower_bound:upper_bound + 1]
            ax.text(x_axis[0], 0.5 * peak_sum / 100. + peak.max(), normalized_integrations[k],
                    fontsize=8)
        if show_plot:
            plt.show()


class MultipleNMRInitializer:
    def __init__(self, file_path):
        self.file_path = file_path  # Path to folder containing multiple NMR files
        self.singlenmrinitializer_dic = 0
        self.to_singlenmeinitializer_dic()
        """"""

    def to_sorted_dict_keys(self, dict):  # to organize the ylabels
        keys = list(dict.keys())
        keys.sort()
        return keys

    def to_singlenmeinitializer_dic(self):
        """Return: dictionary {filename: SingleNMRInitializer object, ... } for all files in directory"""
        directory = os.fsencode(self.file_path)
        directory_str = directory.decode()
        singlenmrinitializer_dic = {}
        for file in os.listdir(directory):
            filename = os.fsdecode(file)
            if filename == ".DS_Store":
                pass
            else:
                pass
                file_path = directory_str + filename
                singlenmrinitializer_dic[filename] = SingleVarian(file_path)
        self.singlenmrinitializer_dic = singlenmrinitializer_dic
        return singlenmrinitializer_dic

    def to_3d_plot(self, show_plot=False, save_fig=False):
        fig = plt.figure(figsize=(10, 6),)
        ax = fig.add_subplot(**{"projection": "3d"})
        singlenmrinitializer_dic = self.singlenmrinitializer_dic
        sorted_keys = self.to_sorted_dict_keys(singlenmrinitializer_dic)

        def get_color_gradient_list(gradient_length, colormap_key):
            grey_colormap = colormaps[colormap_key]
            rgba_list = []
            for item in np.linspace(0, 1, gradient_length*2):
                rgba_list.append(grey_colormap(item))
            return rgba_list

        color_list = get_color_gradient_list(len(sorted_keys), "RdPu")  # change colormap here
        y = len(sorted_keys)-1
        for key, color in zip(sorted_keys[::-1], color_list[len(sorted_keys)-1:]):
            x_values = singlenmrinitializer_dic[key].data_ppm
            y_value = y
            z_values = singlenmrinitializer_dic[key].processed_data
            ax.plot(x_values, y_value, z_values, color=color)
            y -= 1

        def adjust_x_coordinate_text(x_values):
            min_x_values = np.min(x_values)
            step = 1.5  # smaller step will move further from axis
            x_step = (min_x_values + np.min(x_values)) / step
            x_coordinate = min_x_values + x_step
            return x_coordinate

        ax.set_yticks(ticks=np.arange(0, len(sorted_keys), 1), labels=sorted_keys, ha='left')
        ax.tick_params(axis='y', pad=10)
        ax.text(adjust_x_coordinate_text(singlenmrinitializer_dic[sorted_keys[0]].data_ppm),
                0.0, 0.0, "Sample#-example#-hr", zdir='y')  # hard to make ax.set_ylabel look good this works ok
        ax.set_xlabel("σ ppm")
        ax.set_zlabel("Intensity")
        ax.invert_xaxis()
        ax.grid(False)
        fig.align_labels()
        if save_fig:
            plt.savefig(f"{sorted_keys[0]} - {sorted_keys[len(sorted_keys)-1]} to_3d_plot() fig.png", bbox_inches='tight', pad_inches=0.25, dpi=300)
        if show_plot:
            plt.show()
        plt.close()

    def to_integrations_df(self):
        singlenmrinitializer_dic = self.singlenmrinitializer_dic
        sorted_keys = self.to_sorted_dict_keys(singlenmrinitializer_dic)
        integrations_list = []
        data_index_list = []
        experiment_id_list = []
        ppm_list = []
        for key in sorted_keys:
            normalized_integrations = singlenmrinitializer_dic[key].normalized_integrations
            data_ppm_values = singlenmrinitializer_dic[key].data_ppm
            integrations_list += list(normalized_integrations.values())
            data_index_list += list(normalized_integrations.keys())
            ppm_list += list(np.round(data_ppm_values[list(normalized_integrations.keys())], decimals=3))
            experiment_id_list += [key] * len(list(normalized_integrations.values()))
        df = pd.DataFrame({"experiment_id": experiment_id_list, "data_index": data_index_list, "data_ppm": ppm_list, "integration": integrations_list})
        return df

    def to_internal_std_correction(self, set_peak_width_ppm=0.01, integration_value_int_std_peak=1, internal_std_peak_ppm=1.0, peak_to_analyze=0.0):
        df = self.to_integrations_df()
        df.sort_values(by="data_index", inplace=True)

        def make_ppm_groups(set_peak_width_ppm):
            data_indexes = list(df["data_ppm"])
            group = 1
            group_list = []
            for i in range(len(data_indexes)):
                if i <= len(data_indexes)-2:
                    group_list.append(group)
                    if not data_indexes[i+1]+set_peak_width_ppm >= data_indexes[i] >= data_indexes[i+1]-set_peak_width_ppm:
                        group += 1
                else:
                    if not data_indexes[i-1]-set_peak_width_ppm >= data_indexes[i] >= data_indexes[i-1]+set_peak_width_ppm:
                        group += 1
                        group_list.append(group)
                    else:
                        group_list.append(group)
            return group_list
        df["ppm_groups"] = make_ppm_groups(set_peak_width_ppm)

        def get_data_index_ppm(df):
            closest_values_internal_std = abs(df["data_ppm"] - internal_std_peak_ppm)
            internal_std_ppm_group_index = closest_values_internal_std.idxmin()
            internal_std_group = df.loc[internal_std_ppm_group_index, "ppm_groups"]

            closest_values_peak_to_analyze = abs(df["data_ppm"] - peak_to_analyze)
            peak_to_analyze_ppm_group_index = closest_values_peak_to_analyze.idxmin()
            peak_to_analyze_group = df.loc[peak_to_analyze_ppm_group_index, "ppm_groups"]

            return internal_std_group, peak_to_analyze_group

        def to_internal_standard_correction(df, integration_value_int_std_peak):
            internal_std_group, peak_to_analyze_group = get_data_index_ppm(df)  # add func to automate which group is selected based of user input for ppm
            df_copy = df[df["ppm_groups"] == internal_std_group].copy()
            df_copy["internal_std_correction"] = integration_value_int_std_peak / df_copy["integration"]
            df["internal_std_correction"] = np.nan
            for exp_id in list(df_copy["experiment_id"]):
                indexes = df[df["experiment_id"] == exp_id].index
                value_index = df_copy[df_copy["experiment_id"] == exp_id]["internal_std_correction"].index
                value = df_copy._get_value(value_index[0], "internal_std_correction")
                value = [value] * len(list(indexes))
                df.loc[indexes, "internal_std_correction"] = value
            df["integration_internal_std"] = df["integration"] * df["internal_std_correction"]

            return df, internal_std_group, peak_to_analyze_group

        return to_internal_standard_correction(df, integration_value_int_std_peak)

    def to_concentration_over_time(self, set_peak_width_ppm=0.01, integration_value_int_std_peak=1, internal_std_peak_ppm=1.0, peak_to_analyze=0.0, number_of_protons_peak_to_analyze=1):
        df, internal_std_group, peak_to_analyze_group = self.to_internal_std_correction(set_peak_width_ppm, integration_value_int_std_peak, internal_std_peak_ppm, peak_to_analyze)  # data index group needs to be changed to the signal of intrist not the internal std
        df_copy = df[df["ppm_groups"] == peak_to_analyze_group].copy()
        df_copy["integration_mol"] = df_copy["integration_internal_std"] / number_of_protons_peak_to_analyze
        df_copy.sort_values(by="experiment_id", inplace=True)
        return df_copy

    def to_plot_concentration_over_time(self, set_peak_width_ppm=0.01, integration_value_int_std_peak=1, internal_std_peak_ppm=1.0,
                                        peak_to_analyze=0.0, number_of_protons_peak_to_analyze=1, save_fig=False, show_plot=False):
        df = self.to_concentration_over_time(set_peak_width_ppm, integration_value_int_std_peak, internal_std_peak_ppm,
                                             peak_to_analyze, number_of_protons_peak_to_analyze)
        reaction_times = []
        for exp_id in list(df["experiment_id"]):
            time = exp_id.split("-")
            reaction_times.append(float(time[-1]))
        plt.plot(reaction_times, df["integration_mol"], 'ro')
        y_min = 0
        y_max = df["integration_mol"].max()
        x_max = max(reaction_times)
        plt.axis((-10, x_max+10, y_min, y_max+1))
        plt.ylabel("Relative Concentration")
        plt.xlabel("Reaction Time")
        plt.suptitle(f"integrations at {round(df['data_ppm'].mean(), 2)}ppm\n standard at {internal_std_peak_ppm} ppm\n"
                     f"{list(df['experiment_id'])[0]} to {list(df['experiment_id'])[len(list(df['experiment_id']))-1]}", y=1.05)
        if save_fig:
            plt.savefig(f"{list(df['experiment_id'])[0]} to {list(df['experiment_id'])[len(list(df['experiment_id']))-1]}"
                        f" integrations at ppm = {round(df['data_ppm'].mean(), 2)}.png",
                        bbox_inches='tight', pad_inches=0.25, dpi=300)
        if show_plot:
            plt.show()
        plt.close()
