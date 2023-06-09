function [ ] = plot_states( t, sim, species_names, varargin)

experiment = '';
if length(varargin) > 0
    experiment = varargin{1};
end

species_of_interest = {'optoSos', 'optoSos_acc', 'optoSos_star', 'FgfR', 'FgfR_exc', 'FgfR_star', 'FgfR_endo', 'FgfR_auto', 'OptoFgfR_star', ...
    'Sos_star', 'NfbSos_star', 'Ras_star', 'Raf_star', 'NfbRaf_star', 'Nfb_star','Mek_star', ...
    'Erk_star',  'NfbCal_star',  'KTR_star', 'Ktr_star',  'Iff_star', 'Chan_acc', ...
    'ChanOut_acc', 'ChanIn_acc', 'Cal', 'Cal_star', 'Calcin_star'};

title_index = 1;
for i = 1: length(species_of_interest)
    for j = 1: length(species_names)
        if strcmp(species_names{j}, species_of_interest{i})
            string = strsplit(species_of_interest{i}, '_');
            title_tmp =['$\textnormal{', string{1}, '}$'];
            if length(string) > 1
            if strcmp(string{2}, 'star')
                title_tmp = ['$\textnormal{', string{1}, '}^*$'];
            else
                title_tmp =  ['$\textnormal{', string{1}, '}_{\textnormal{', string{2} '}}$'];
            end
            end
            titles{title_index} =title_tmp;
            species_indices(title_index) = j;
            title_index = title_index+1;
        end
    end
end

num_plots = length(titles);
num_cols = 2;
num_species = length(species_names);

for i = 1 : num_plots
    subplot(ceil(num_plots /num_cols), num_cols, i);
    hold on;
    plot(t, sim(species_indices(i) : num_species : end , :), 'LineWidth', 2);
    title([titles{i}, ' ', experiment], 'Interpreter', 'Latex');
end
end

