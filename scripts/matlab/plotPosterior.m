function [] = plotPosterior(summary_file_name, varargin)
%% The input is
%   summary_file_name: The path and name of the model summary file. The
%                       path can be absolute or relative
%                       (./results/birth_death_model_summary.txt for
%                       instance. The model summary file must be in the
%                       same folder as all other results files

%% Optional input: 
%   plot_scatter:   true or false, indicates if a scatterplot of pairwise
%                   posterior is to be plotted. 
%   plot_level_set: true or false, indicates if the last level set instead
%                   of the posterior should be plotted
%   plot_highest_corr: int, plots only the indicated number of highest
%                   correlated parameters

selected_params = {};
ignored_parameters = {'ktr_init', 'kin_init', 'phos_init'};


clear global
max_it = 5000000;
folder_index = strfind(summary_file_name, '/');
file_name_index = strfind(summary_file_name, '_model_summary.txt');
plot_scatter = false;
plot_level_set = false;
plot_highest_corr = -1;

if(length(varargin) > 0)
    for i = 1: 2: length(varargin)
        option_name = varargin{i};
        if strcmp(option_name, 'plot_scatter')
            plot_scatter = varargin{i+1};
        elseif strcmp(option_name, 'plot_level_set')
            plot_level_set = varargin{i+1};
        elseif strcmp(option_name, 'plot_highest_corr')
            plot_highest_corr = varargin{i+1};
            
        end
    end
end
if isempty(file_name_index)
    fprintf('The provided file must be the *_model_summary.txt file. It must be in the same folder as the other output files ( *_times.txt, _latent_states*.txt, *_measurements.txt)' );
    return;
end

folder_name = summary_file_name(1:folder_index(end));
file_name = [folder_name, summary_file_name(folder_index(end) +1 : file_name_index)];


[param_names, species_names, scales, bounds] = readModelDescription(summary_file_name);


if ~plot_level_set
    [posterior, weights] = getPosteriorParticles(summary_file_name, max_it);
else
    results_full_name_split =  strsplit(summary_file_name, '_model_summary');
    results_full_name = results_full_name_split{1};
    results_file_name_split = strsplit(results_full_name, '/');
    results_file_name = results_file_name_split{end};
    results_folder_name = results_full_name(1: end - length(results_file_name));
    
    
    results_files = dir(strcat(results_folder_name));
    
    file_like = regexpi({results_files.name}, strcat(results_file_name, '_live_points_[0-9]*_log_likelihoods.txt'),'match');
    
    file_like = [file_like{:}];max_nbr_like = -1;
    max_ind_like = -1;
    for i = 1 : length(file_like)
        str = strsplit(file_like{i}, '_');
        nbr = str2num(str{4});
        if(nbr > max_nbr_like)
            max_nbr_like = nbr; max_ind_like = i;
        end
    end
    posterior = dlmread([results_folder_name, '/', results_file_name, '_live_points_', num2str(max_nbr_like), '.txt']);
    weights = ones(size(posterior, 1), 1);
    
end

if ~isempty(ignored_parameters)
    selected_params = param_names;
    indices = [];
    for j= 1: length(ignored_parameters)
        index = find(contains(param_names,ignored_parameters{j}));
        indices = [indices, index];
    end
    selected_params(indices) = [];
end

if plot_highest_corr > 0
    corr_mat =  corr(posterior);
    corr_mat = tril(corr_mat, -1);
    corr_mat = abs(corr_mat);
    max_elements = maxk(corr_mat(:), plot_highest_corr);
    max_indices = [];
    for i = 1: plot_highest_corr
        index = find(corr_mat(:) == max_elements(i));
        
        s = size(corr_mat);
        [row, col] = ind2sub(s,index);
        max_indices = [max_indices; row, col];
        selected_params = [selected_params, param_names{max_indices(i, 1)}];
        selected_params = [selected_params, param_names{max_indices(i, 2)}];
    end
    selected_params = unique(selected_params);
end

if ~isempty(selected_params)
    param_indices = [];
    for i = 1: length(param_names)
        for j = 1 : length(selected_params)
            if strcmp(param_names{i}, selected_params{j})
                param_indices = [param_indices, i];
            end
        end
    end
    posterior = posterior(:, param_indices);
    param_names = param_names(param_indices);
    scales = scales(param_indices);
    bounds = bounds(param_indices, :);
end

for i = 1: length(scales)
    if strcmp( scales{i}, 'log')
        posterior(:, i) = log10(posterior(:, i));
    end
end

num_params = length(param_names);
num_cols = 1;
if num_params > 1
    num_cols = 3;
end
num_rows = ceil(num_params / 3);

samples = randsample(1:size(posterior, 1), 100000, true, weights);

figure();
prior_interval = bounds;
for i = 1 : num_params
    subplot(num_rows, num_cols, i);
    [f, xi] =   ksdensity(posterior(:, i), 'weights', weights);%, 'Bandwidth', 0.5);
    plot(xi, f);
    hold on;
    [counts, bin] = hist(posterior(samples, i));
    %     bar(bin, (counts / (max(counts))) * max(f), 'FaceColor', [0, .45, .74]);
    bar(bin, (counts / (max(counts))) * max(f), 'FaceColor', [0.85, .33, .1]);
    %     bar(bin, (counts / (max(counts))) * max(f), 'FaceColor', [0.93, .69, .13]);
    hold on;
    if strcmp(scales{i}, 'log') == 1
        title(['log(', param_names{i}, ')']);
        xlim(log10(bounds(i, :)));
        prior_interval(i, : ) = log10(bounds(i, :));
    else
        title(param_names{i});
        xlim(bounds(i, :));
    end
end

%% plot scatter plot
if plot_scatter
    post = posterior(samples, :)';
    figure();
    for i = 1: num_params
        for j = 1: i
            if(i == j)
                subplot(num_params, num_params, (i-1)*num_params + j);
                ksdensity(post(i, :));
                %             xlabel(param_names{i}, 'Interpreter', 'Latex');
                xlim(prior_interval(i, :))
            else
                subplot(num_params, num_params, (i-1)*num_params + j);
                grid = 20;%256;   %refinement of map
                points = [post(j, :)', post(i, :)'];
                minvals = [prior_interval(i, 1), prior_interval(j, 1)];%min(points);
                maxvals = [prior_interval(i, 2), prior_interval(j, 2)];%max(points);
                rangevals = maxvals - minvals;
                xidx = 1 + round((points(:,1) - minvals(2)) ./ rangevals(2) * (grid-1));
                yidx = 1 + round((points(:,2) - minvals(1)) ./ rangevals(1) * (grid-1));
                density = accumarray([yidx, xidx], 1, [grid,grid]);  %note y is rows, x is cols
                imagesc(density, 'xdata', [minvals(2), maxvals(2)], 'ydata', [minvals(1), maxvals(1)]);
                %             imagesc(density, 'xdata', [-1 1], 'ydata', [-5 5]);
                set(gca,'YDir','normal')
                
                if j == 1
                    ylabel(param_names{i}, 'Interpreter', 'Latex');
                end
            end
            if i == num_params
                xlabel(param_names{j}, 'Interpreter', 'Latex');
            end
        end
    end
end
end