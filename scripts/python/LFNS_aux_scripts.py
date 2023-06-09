import re
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import pandas as pd
import csv
import os
from matplotlib.ticker import FuncFormatter



def read_model_description(path_to_model_description):
    if not os.path.isfile(path_to_model_description):
        print("Error: The provided model description file cannot be found.")
        return

    with open(path_to_model_description, 'r') as file:
        model_summary = file.read()

    # Split the model summary into the relevant sections
    model_sections = re.split(
        r'[-]{3,}\s+Model Settings\s+[-]{3,}|[-]{3,}\s+Model\s+[-]{3,}|[-]{3,}\s+Initial Values\s+[-]{3,}',
        model_summary)
    lfns_info_substr = model_sections[0]
    model_setting_substr = model_sections[1]
    model_info_substr = model_sections[2]
    initial_value_substr = model_sections[3]

    # Extract experiments if they exist
    experiments_line = re.search(r'Experiments for.*:(.*?)(?=\n|$)', lfns_info_substr)
    if experiments_line:
        experiments = [x.strip() for x in experiments_line.group(1).split(',') if x.strip()]

    # Extract parameter bounds
    bounds_info = re.findall(r'Name\s+Bounds\s+Scale\n(.*)', lfns_info_substr, re.DOTALL)
    param_names, scales, bounds = [], [], []
    if bounds_info:
        bounds_lines = [x for x in bounds_info[0].split('\n') if x]
        for line in bounds_lines:
            entries = line.strip().split()
            if len(entries) > 2:
                param_names.append(entries[0])
                scales.append(entries[-1])
                bounds.append([float(entries[-3][1:-1]), float(entries[-2][:-1])])
    else:
        param_info = re.findall(r'Model parameters\s\n(.*)', model_setting_substr, re.DOTALL)
        param_lines = [x for x in param_info[0].split('\n') if x]
        for line in param_lines:
            entries = line.strip().split()
            if len(entries) == 1:
                param_names.append(entries[0])
    # Extract species names
    species_names = re.findall(r'Species:(.*?)(?=\n|$)', model_info_substr)
    if species_names:
        species_names = species_names[0].strip().split()

    # Extract provided parameters and file if they exist
    provided_params, provided_params_file = [], None
    if 'provided the parameters' in model_summary:
        provided_params_info = re.findall(r'provided the parameters(.*?)(?=Parameters to be sampled:|$)', model_summary,
                                          re.DOTALL)
        if provided_params_info:
            provided_params = provided_params_info[0].strip().split()
            provided_params_file = \
            re.findall(r'will be read from the file(.*?)(?=Parameters to be sampled:|$)', model_summary, re.DOTALL)[
                0].strip()

    model_summaries = {"param_names": param_names, "species_names": species_names, "scales": scales,
                     "bounds": bounds, "experiments": experiments, "provided_params": provided_params,
                     "provided_params_file": provided_params_file }
    return model_summaries


def plot_posterior(summary_file_name, plot_scatter=False, plot_highest_corr=-1, plot_level_set = False,
                  ignored_parameters=[]):
    max_it = 5000000

    model_summary = read_model_description(summary_file_name)

    if plot_level_set:
        folder_name = Path(summary_file_name).parent
        file_name = Path(summary_file_name).name.replace('_model_summary.txt', '')
        logs = pd.read_csv(folder_name / f"{file_name}_log_file.txt", sep='\t', skiprows=1, header=None)
        max_iteration_nbr = logs.shape[0]
        posterior = np.loadtxt(folder_name / f"{file_name}_live_points_{max_iteration_nbr}.txt")
        weights = np.ones(posterior.shape[0])
    else:
        posterior, weights = get_posterior_particles(summary_file_name, max_it)
    selected_params = [p for p in model_summary["param_names"] if p not in ignored_parameters]

    if selected_params:
        param_indices = [i for i, p in enumerate(model_summary["param_names"]) if p in selected_params]
        posterior = posterior[:, param_indices]
        model_summary["param_names"] = [model_summary["param_names"][i] for i in param_indices]
        model_summary["scales"] = [model_summary["scales"][i] for i in param_indices]
        model_summary["bounds"] = [model_summary["bounds"][i] for i in param_indices]


    for i in range(len(model_summary["scales"])):
        if model_summary["scales"][i] == 'log':
            model_summary["bounds"][i]  = np.log10(model_summary["bounds"][i])
            posterior[:, i] = np.log10(posterior[:, i])

    num_params = len(model_summary["param_names"])
    num_cols = 1 if num_params <= 1 else 3
    num_rows = int(np.ceil(num_params / 3))

    weights_norm = weights / sum(weights)
    samples = np.random.choice(posterior.shape[0], 1000, p=weights_norm)
    posterior_samples = posterior[samples, ]

    fig, axs = plt.subplots(num_rows, num_cols)
    fig.tight_layout()
    axs = axs.ravel()  # Flatten the array of subplots

    for i in range(num_params):
        ax = axs[i]
        non_zero_indices = weights_norm > 0
        # print(posterior[non_zero_indices, i])
        # print(weights_norm[non_zero_indices])
        density = gaussian_kde(posterior[non_zero_indices, i], weights= weights_norm[non_zero_indices])
        xs = np.linspace(*model_summary["bounds"][i], 200)
        ax.plot(xs, density(xs))
        counts, bins = np.histogram(posterior_samples[:, i])
        counts = counts * (max(density(xs)) / max(counts))
        ax.stairs(counts, bins)
        ax.set_title('log({})'.format(model_summary["param_names"][i]) if model_summary["scales"][i] == 'log' else model_summary["param_names"][i])

    if plot_scatter:
        pass

    plt.show()


def get_posterior_particles(summary_file_name, max_it=1e6):
    if not os.path.isfile(summary_file_name):
        print("Error: The provided summary file cannot be found.")
        return
    folder_name = Path(summary_file_name).parent
    file_name = Path(summary_file_name).name.replace('_model_summary.txt', '')

    # Read in the data
    posterior_particles = np.loadtxt(folder_name / f"{file_name}_posterior.txt")
    log_weights = np.loadtxt(folder_name / f"{file_name}_posterior_log_weights.txt")

    # Compute the weights
    full_weights = np.exp(log_weights - np.max(log_weights))

    return posterior_particles, full_weights


def plot_system(summary_file_name, data_file_names = [], time_file_names = []):
    folder_index = summary_file_name.rfind('/')
    file_name_index = summary_file_name.rfind('_model_summary.txt')

    measure_datas = []
    time_datas = []
    if len(data_file_names)>0:
        for data_file_name in data_file_names:
            with open(data_file_name, 'r') as file:
                csvreader = csv.reader(file)
                for row in csvreader:
                    measure_data = row
                    measure_data = list(map(float, measure_data))
            measure_datas.append(measure_data)
        for time_file_name in time_file_names:
            with open(time_file_name, 'r') as file:
                csvreader = csv.reader(file)
                for row in csvreader:
                    time_data = row
                    time_data = list(map(float, time_data))
            time_datas.append(time_data)

    if file_name_index == -1:
        print('The provided file must be the *_model summary.txt file. It must be in the same folder as the other output files ( *_times.txt, _latent_states*.txt, *_measurements.txt)' )
        return

    if folder_index != -1:
        folder_name = summary_file_name[:folder_index+1]
        file_name = folder_name + summary_file_name[folder_index+1:file_name_index]
    else:
        folder_name = './'
        file_name = folder_name + summary_file_name[:file_name_index]

    model_summary = read_model_description(summary_file_name)
    for count, experiment in enumerate(model_summary["experiments"]):
        latent_states_file = file_name + "_" + experiment + '_latent_states.txt'
        measurement_states_file = file_name + "_" + experiment + '_measurements.txt'
        times_file = file_name + '_times.txt'
        t = np.loadtxt(times_file, delimiter=',')
        latent_states = np.loadtxt(latent_states_file, delimiter=',')
        measurement = np.loadtxt(measurement_states_file, delimiter=',')

        num_states = len(model_summary["species_names"])
        num_simulations = latent_states.shape[0] // num_states

        num_cols = 1
        if num_states > 1:
            num_cols = 2
        num_rows = int(np.ceil(num_states / 2))

        if num_simulations == 1:
            plot_fig_ofset =((count-1)* 2 )
        else:
            plot_fig_ofset =((count-1)* 3 )

        plt.figure(plot_fig_ofset + 1,  figsize=(6, 3))
        if num_simulations == 1:
            plt.plot(t, measurement, '--')
            if len(time_datas) > 0:
                plt.plot(np.array(time_datas[count]), np.array(measure_datas[count]),color="red")
            plt.xlabel('time')
            plt.ylabel('measurement')
            plt.title('Measurements ' + experiment)
        else:
            plt.subplot(1, 2, 1)
            [plt.plot(t, measurement[row,], 'o', markersize=0.5)for row in range(num_simulations)]
            if len(time_datas) > 0:
                plt.plot(np.array(time_datas[count]), np.array(measure_datas[count]),color="red")
            plt.xlabel('time')
            plt.ylabel('measurement')
            plt.title('Measurements ' + experiment)

            plt.subplot(1, 2, 2)
            plt.plot(t, np.mean(measurement, axis=0), 'o', markersize=0.5)
            if len(time_datas) > 0:
                if np.array(measure_datas[count]).ndim>1:
                    plt.plot(np.array(time_datas[count]), np.mean(np.array(measure_datas[count]), axis=0),color="red")
                else:
                    plt.plot(np.array(time_datas[count]), np.array(measure_datas[count]),color="red")

            plt.xlabel('time')
            plt.ylabel('measurement')
            plt.title('Mean measurements ' + experiment)
        plt.tight_layout()

        plt.figure(plot_fig_ofset + 2, figsize=(3*num_cols, 3*num_rows))
        for i in range(num_states):
            plt.subplot(num_rows, num_cols, i+1)
            [plt.plot(t, latent_states[(num_states *row) + i::num_states *row+1 , :][0], linewidth=2) for row in range(num_simulations) ]
            plt.xlabel('time')
            plt.ylabel('state')
            plt.title(model_summary["species_names"][i] + ' ' + experiment)
        plt.tight_layout()

        if num_simulations > 1:
            plt.figure(plot_fig_ofset + 3, figsize=(6, 3))
            for i in range(num_states):
                plt.subplot(num_rows, num_cols, i+1)
                plt.plot(t, np.mean(latent_states[i::num_states, :], axis=0), linewidth=2)
                plt.xlabel('time')
                plt.ylabel('state')
                plt.title('mean ' + model_summary["species_names"][i] + ' ' + experiment)
            plt.tight_layout()
    plt.show()


def plotBE(log_file_name):
    logs = pd.read_csv(log_file_name, sep='\t', skiprows=1, header=None)
    max_iteration_nbr = logs.shape[0]

    log_zd = logs.iloc[:, 6]
    log_zl = logs.iloc[:, 8]
    log_ztot = logs.iloc[:, 4]

    log_var_zd = logs.iloc[:, 7]
    log_var_zl = logs.iloc[:, 9]
    log_var_ztot = logs.iloc[:, 5]
    log_var_min = logs.iloc[:, 10]
    indices = np.arange(max_iteration_nbr)

    log_std_zd = 0.5 * log_var_zd
    max_log = np.maximum(log_zd, log_std_zd)
    zd_down = np.log(np.maximum(np.exp(log_zd - max_log) - np.exp(log_std_zd - max_log), 0)) + max_log
    zd_up = np.log(np.exp(log_zd - max_log) + np.exp(log_std_zd - max_log)) + max_log
    zd_down = np.where(np.isinf(zd_down), np.min(log_zd), zd_down)

    plt.fill_between(indices, zd_down, zd_up, color=[.93, .84, .84])

    log_std_zl = 0.5 * log_var_zl
    max_log = np.maximum(log_zl, log_std_zl)
    zl_down = np.log(np.maximum(np.exp(log_zl - max_log) - np.exp(log_std_zl - max_log), 0)) + max_log
    zl_up = np.log(np.exp(log_zl - max_log) + np.exp(log_std_zl - max_log)) + max_log
    zl_down = np.where(np.isinf(zl_down), np.min(log_zd), zl_down)

    plt.fill_between(indices, zl_down, zl_up, color=[.76, .87, .78])

    log_std_ztot = 0.5 * log_var_ztot
    max_log = np.maximum(log_ztot, log_std_ztot)
    zLFNS_down = np.log(np.maximum(np.exp(log_ztot - max_log) - np.exp(log_std_ztot - max_log), 0)) + max_log
    zLFNS_up = np.log(np.exp(log_ztot - max_log) + np.exp(log_std_ztot - max_log)) + max_log
    zLFNS_down = np.where(np.isinf(zLFNS_down), np.min(log_zd), zLFNS_down)

    plt.fill_between(indices, zLFNS_down, zLFNS_up, color=[.73, .83, .96])
    plt.plot(log_zd, linewidth=1)
    plt.plot(log_zl, linewidth=1)
    plt.plot(log_ztot, linewidth=1)
    plt.show()

    plt.figure()

    lfns_error = logs.iloc[:, -2]
    max_error = logs.iloc[:, -3]
    times = logs.iloc[:, 3]
    acceptance = logs.iloc[:, 2]

    plt.figure(figsize=(20, 5))
    plt.subplot(1, 4, 1)
    plt.plot(indices, max_error, linewidth=2)
    plt.plot(indices, lfns_error, linewidth=2)
    plt.title('Error estimate')
    plt.legend(['$\log(\Delta_{max}^m)$', '$\log(\Delta_{LFNS}^m)$'], loc='best')
    plt.xlabel('Iteration Nbr m')
    plt.ylabel('$\log(\Delta^m)$')
    plt.grid(True)
    plt.xlim(indices[0], indices[-1] + 1)

    plt.subplot(1, 4, 2)
    plt.plot(indices, 0.5 * log_var_ztot, linewidth=2, color=[.73, .83, .96])
    plt.plot(indices, 0.5 * log_var_min, linewidth=2, color=[.31, .4, .58])
    plt.title('Variance estimate')
    plt.legend(['$\log(\hat{\sigma}_{tot}^{2m})$', '$\log(\hat{\sigma}_{min}^{2m})$'], loc='best')
    plt.xlabel('Iteration Nbr m')
    plt.ylabel('$\log(\sigma^{2m})$')
    plt.grid(True)
    plt.xlim(indices[0], indices[-1] + 1)

    plt.subplot(1, 4, 3)
    norm_const = log_ztot.iloc[-1] if np.exp(log_ztot.iloc[-1]) > 0 else 0

    # Calculate z and std normalized values
    zd_norm = np.exp(log_zd - norm_const)
    log_var_zd_norm = log_var_zd - 2 * norm_const
    std_zd_norm = np.exp(0.5 * log_var_zd_norm)

    # Calculate down and up values
    zd_down = np.maximum(zd_norm - std_zd_norm, 0)
    zd_up = zd_norm + std_zd_norm

    zl_norm = np.exp(log_zl - norm_const)
    log_var_zl_norm = log_var_zl - 2 * norm_const
    std_zl_norm = np.exp(0.5 * log_var_zl_norm)
    zl_down = np.maximum(zl_norm - std_zl_norm, 0)
    zl_up = zl_norm + std_zl_norm

    ztot_norm = np.exp(log_ztot - norm_const)
    log_var_ztot_norm = log_var_ztot - 2 * norm_const
    std_ztot_norm = np.exp(0.5 * log_var_ztot_norm)
    zLFNS_down = np.maximum(ztot_norm - std_ztot_norm, 0)
    zLFNS_up = ztot_norm + std_ztot_norm

    log_var_min_norm = log_var_min - 2 * norm_const
    std_std_min_norm = np.exp(0.5 * log_var_min_norm)
    zmin_down = np.maximum(ztot_norm - std_std_min_norm, 0)
    zmin_up = ztot_norm + std_std_min_norm

    plt.fill_between(indices, zd_down, zd_up, color=[.93, .84, .84])
    plt.fill_between(indices, zl_down, zl_up, color=[.76, .87, .78])
    plt.fill_between(indices, zLFNS_down, zLFNS_up, color=[.73, .83, .96])
    plt.fill_between(indices, zmin_down, zmin_up, color=[.31, .4, .58])

    plt.plot(zd_norm, linewidth=2, color=[.49, .18, .56])
    plt.plot(zl_norm, linewidth=2, color=[.47, .67, .19])
    plt.plot(ztot_norm, linewidth=2, color=[0.3, .75, .93])

    plt.title('BE estimate')
    plt.legend(['$\pm\hat{\sigma}_\mathcal{D}^{m}$', '$\pm\hat{\sigma}_\mathcal{L}^{m}$',
                '$\pm\hat{\sigma}_{tot}^{m}$', '$\pm\hat{\sigma}_{min}^{m}$',
                 '$\hat{Z}_{\mathcal{D}}$', '$\hat{Z}_{\mathcal{L}}$', '$\hat{Z}_{tot}$',],
               prop={"size": 10}, loc='best')
    plt.ylabel('$\hat{Z}$')
    plt.xlabel('Iteration Nbr m')
    plt.xlim([indices[0], indices[-1] + 1])

    plt.subplot(1, 4, 4)
    plt.grid(True)
    plt.plot(indices, acceptance, linewidth=2)
    plt.ylabel('Acceptance rate')
    plt.xlabel('Iteration Nbr m')


    ax2 = plt.gca().twinx()  # Create a second y-axis
    cum_times = np.cumsum(times)
    if cum_times.iloc[-1] > 2 * 3600:
        cum_times = cum_times / 3600
        ax2.set_ylabel('Cumulative time in hours')
    elif cum_times.iloc[-1] > 2 * 60:
        cum_times = cum_times / 60
        ax2.set_ylabel('Cumulative time in minutes')
    else:
        ax2.plot(indices, cum_times, linewidth=2, color = "#ff7f0e" )
        ax2.set_ylabel('Cumulative time in seconds')

    plt.title('Computational effort')
    plt.show()


