function [param_names, species_names, scales, bounds, experiments, provided_params, provided_params_file] = readModelDescription(path_to_model_description)

model_summary = fileread(path_to_model_description);
model_key = '-----------   Model   -----------';
initial_value_key = '----------   Initial Values   ---------- ';

lfns_substring_index = strfind(model_summary,model_key);
lfns_info_substr = model_summary(1 : lfns_substring_index);

model_substring_index = strfind(model_summary,initial_value_key);
model_info_substr = model_summary(lfns_substring_index +1 : model_substring_index);

model_conf_key = 'Experiments for';
model_configurations_index = strfind( lfns_info_substr, model_conf_key);
experiments = '';
if(length(model_configurations_index) > 0)
    model_conf_lines = splitlines(lfns_info_substr(model_configurations_index : end));
    model_conf_lines = model_conf_lines{1};
    model_conf_lines = strsplit(model_conf_lines, ':');
    model_conf_lines = model_conf_lines{2};
    experiments = strsplit(model_conf_lines, ', ');
    for i = 1: length(experiments)
        experiments{i} = strip(experiments{i} );
        if( length(experiments{i}) == 0)
            experiments(i) = [];
        end
    end
end

bounds_key = 'Name                  Bounds  Scale';

bounds_index = strfind( lfns_info_substr, bounds_key);
bounds_lines = splitlines(lfns_info_substr(bounds_index : end));
bounds_lines(1) = [];
for i = 1 : length(bounds_lines)
    if(length(bounds_lines{i}) == 0)
        bounds_lines(i:end) = [];
        break;
    end
end


param_names = {};
scales= {};
bounds = [];
for line_nbr = 1 : length(bounds_lines)
    entries = strsplit(strip(bounds_lines{line_nbr}), ' ');
    
    if(length(entries) >2)
        param_name = entries{1};
        scale = entries{end};
        r_bound = entries{end-1};
        r_bound = str2num(r_bound(1:end-1));
        
        l_bound =  entries{end-2};
        l_bound = str2num(l_bound(2:end-1));
        bound = [l_bound, r_bound];
        
        param_names = [param_names, param_name];
        scales = [scales; scale];
        bounds = [bounds; bound];
    end
end

species_key = 'Species:';
species_index = strfind( model_info_substr, species_key);
species_lines = splitlines(model_info_substr(species_index : end - 1));
species_line = species_lines{1};

species_names = strsplit(strip(species_line), ' ');
species_names(1) = [];


provided_params_key = 'provided the parameters';
provided_params_index = strfind(model_summary, provided_params_key);

provided_params_file = [];
provided_params = [];
if ~isempty(provided_params_index)
    provided_param_file_key = 'will be read from the file';
    provided_param_file_index = strfind(model_summary, provided_param_file_key);
    provided_params_str = model_summary(provided_params_index +24 : provided_param_file_index -3);
    provided_params = strsplit(provided_params_str, ' ');
    
    sampled_params_key= 'Parameters to be sampled:';
    sampled_params_index = strfind(model_summary, sampled_params_key);
    
    provided_params_file = model_summary(provided_param_file_index +27 : sampled_params_index - 2);
    
    
end
end
