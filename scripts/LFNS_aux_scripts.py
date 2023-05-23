import re
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import pandas as pd


def read_model_description(path_to_model_description):
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
        param_info = re.findall(r'Name\s+simulation value\n(.*)', lfns_info_substr, re.DOTALL)
        param_lines = [x for x in param_info[0].split('\n') if x]
        for line in param_lines:
            entries = line.strip().split()
            if len(entries) == 2:
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

    return param_names, species_names, scales, bounds, experiments, provided_params, provided_params_file


def plot_posterior(summary_file_name, plot_scatter=False, plot_highest_corr=-1,
                  ignored_parameters=[]):
    selected_params = []
    max_it = 5000000
    folder_name = Path(summary_file_name).parent
    file_name = Path(summary_file_name).name.replace('_model_summary.txt', '')

    param_names, species_names, scales, bounds, experiments, provided_params, provided_params_file = read_model_description(summary_file_name)

    posterior, weights = get_posterior_particles(summary_file_name, max_it)
    selected_params = [p for p in param_names if p not in ignored_parameters]
    print(selected_params)
    print(bounds)

    if selected_params:
        param_indices = [i for i, p in enumerate(param_names) if p in selected_params]
        posterior = posterior[:, param_indices]
        param_names = [param_names[i] for i in param_indices]
        scales = [scales[i] for i in param_indices]
        bounds = [bounds[i] for i in param_indices]


    for i in range(len(scales)):
        if scales[i] == 'log':
            bounds[i]  = np.log10(bounds[i])
            posterior[:, i] = np.log10(posterior[:, i])

    num_params = len(param_names)
    num_cols = 1 if num_params <= 1 else 3
    num_rows = int(np.ceil(num_params / 3))

    samples = np.random.choice(posterior.shape[0], 100000, p=weights/sum(weights))

    fig, axs = plt.subplots(num_rows, num_cols)
    axs = axs.ravel()  # Flatten the array of subplots

    for i in range(num_params):
        ax = axs[i]
        density = gaussian_kde(posterior[:, i], weights=weights)
        xs = np.linspace(*bounds[i], 200)
        ax.plot(xs, density(xs))
        ax.set_title('log({})'.format(param_names[i]) if scales[i] == 'log' else param_names[i])

    if plot_scatter:
        pass

    plt.show()


def get_posterior_particles(summary_file_name, max_it=1e6):
    folder_name = Path(summary_file_name).parent
    file_name = Path(summary_file_name).name.replace('_model_summary.txt', '')
    print(folder_name)
    print(file_name)

    # Read in the data
    posterior_particles = np.loadtxt(folder_name / f"{file_name}_posterior.txt")
    log_weights = np.loadtxt(folder_name / f"{file_name}_posterior_log_weights.txt")

    # Compute the weights
    full_weights = np.exp(log_weights - np.max(log_weights))

    return posterior_particles, full_weights


import numpy as np
import matplotlib.pyplot as plt

def plot_system(summary_file_name):
    folder_index = summary_file_name.rfind('/')
    file_name_index = summary_file_name.rfind('_model_summary.txt')

    if file_name_index == -1:
        print('The provided file must be the *_model summary.txt file. It must be in the same folder as the other output files ( *_times.txt, _latent_states*.txt, *_measurements.txt)' )
        return

    if folder_index != -1:
        folder_name = summary_file_name[:folder_index+1]
        file_name = folder_name + summary_file_name[folder_index+1:file_name_index]
    else:
        folder_name = './'
        file_name = folder_name + summary_file_name[:file_name_index]

    param_names, species_names, scales, bounds, experiments, provided_params, provided_params_file = read_model_description(summary_file_name)

    for experiment in experiments:
        latent_states_file = file_name + "_" + experiment + '_latent_states.txt'
        measurement_states_file = file_name + "_" + experiment + '_measurements.txt'
        times_file = file_name + '_times.txt'

        t = np.loadtxt(times_file, delimiter=',')
        latent_states = np.loadtxt(latent_states_file, delimiter=',')
        measurement = np.loadtxt(measurement_states_file, delimiter=',')

        num_states = len(species_names)
        num_simulations = latent_states.shape[0] // num_states

        num_cols = 1
        if num_states > 1:
            num_cols = 2
        num_rows = int(np.ceil(num_states / 2))

        plt.figure(1)
        if num_simulations == 1:
            plt.plot(t, measurement, '--')
            plt.xlabel('time')
            plt.ylabel('measurement')
            plt.title('Measurements ' + experiment)
        else:
            plt.subplot(1, 2, 1)
            plt.plot(t, measurement, 'o')
            plt.xlabel('time')
            plt.ylabel('measurement')
            plt.title('Measurements ' + experiment)

            plt.subplot(1, 2, 2)
            plt.plot(t, np.mean(measurement, axis=0), 'o')
            plt.xlabel('time')
            plt.ylabel('measurement')
            plt.title('Mean measurements ' + experiment)

        plt.figure(2)
        for i in range(num_states):
            plt.subplot(num_rows, num_cols, i+1)
            plt.plot(t, latent_states[i::num_states, :][0], linewidth=2)
            plt.xlabel('time')
            plt.ylabel('state')
            plt.title(species_names[i] + ' ' + experiment)

        if num_simulations > 1:
            plt.figure(3)
            for i in range(num_states):
                plt.subplot(num_rows, num_cols, i+1)
                plt.plot(t, np.mean(latent_states[i::num_states, :], axis=0), linewidth=2)
                plt.xlabel('time')
                plt.ylabel('state')
                plt.title('mean ' + species_names[i] + ' ' + experiment)

    plt.show()
