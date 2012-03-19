clear all;
close all;
granule_cell_filenames = {...
    'DATA2/CD_ipsp/*.mat',...
    'DATA2/CDpause_GCs/*.mat',...
    'DATA2/med/*.mat',...
    'DATA2/PCA/*.mat',...
    'DATA2/PCA_CDpause_GCs/*.mat',...
    'DATA2/PCA_late/*.mat',...
    'DATA2/PCA_med/*.mat',...
    'DATA2/PCA_med_late/*.mat'};
%     'dataformodel/CDpause_GCs/*.mat',...
%     'dataformodel/med_gcs/*.mat',...
%     'dataformodel/PCA_CDpause_GCs/*.mat',...
%     'dataformodel/PCA_GCs/*.mat',...
%     'dataformodel/pca_med_gcs/*.mat',...
%     'dataformodel/PCA_w_late/*.mat'};

plasticity_filenames = {...
    'dataformodel/plasticity/20111213/*.mat',...
    'dataformodel/plasticity/20111220/*.mat'};

% make list of files
gc_files = get_file_list(granule_cell_filenames);

plasticity_files = get_file_list(plasticity_filenames);

% find the maximum and minimum times
min_time = -Inf;
max_time = Inf;
for i=1:length(gc_files);
    gc_info = importdata(char(gc_files(i)));
    times = gc_info.start + gc_info.interval*((1:gc_info.length)-1);
    
    min_time = max(min_time, times(1));
    max_time = min(max_time, times(length(times)));
end
for i=1:length(plasticity_files)
    p_info = importdata(char(plasticity_files(i)));
    times = p_info.start + p_info.interval*((1:p_info.length)-1);
    min_time = max(min_time, times(1));
    max_time = min(max_time, times(length(times)));
end
interval = 1e-3;
times = min_time:interval:max_time;
%% create matrices for granule cell voltages and firing rates
% create array of the granule cell voltages
gc_voltages = zeros(length(times), length(gc_files));
for i=1:length(gc_files)
    gc_info = importdata(char(gc_files(i)));
    gc_times = gc_info.start + gc_info.interval*((1:gc_info.length)-1);
    gc_voltages(:,i) = interp1(gc_times, gc_info.values, times);
end

% convert to rates
gc_rates = zeros(size(gc_voltages));
%rate_parameters = struct('r0', 10, 'DeltaV', 5);
rate_parameters = struct('power',1.5,'v0',0);

figure;
for i=1:size(gc_voltages,2)
    baseline = mean(gc_voltages(1:100,i));
    gc_rates(:,i) = convert_voltages_to_spike_rates(...
        gc_voltages(:,i) - baseline, 'POWER', rate_parameters);
    plot(times, gc_rates(:,i));
    hold all;
end
hold off;

%% fit each plasticity trace will least squares
figure;
pinv_gc_rates = pinv(gc_rates);
for i=1:length(plasticity_files)
    p_info = importdata(char(plasticity_files(i)));
    p_times = p_info.start + p_info.interval*((1:p_info.length)-1);
    
    voltage_change = interp1(p_times, p_info.values, times)';

    weight_changes = pinv_gc_rates * voltage_change;
    fitted_voltage_change =  gc_rates * weight_changes;
    subplot(ceil(length(plasticity_files)/3), 3, i);
    plot(times, voltage_change, times, fitted_voltage_change);
    title(p_info.delay);
end

%% fit the learning rule by least squares
learning_window_times = -0.1:0.005:0.1;

% construct vector by concatenating voltage changes at each delay
voltage_changes = [];
delays = [];
for i=1:length(plasticity_files)
    p_info = importdata(char(plasticity_files(i)));
    p_times = p_info.start + p_info.interval*((1:p_info.length)-1);
    voltage_changes = ...
        [voltage_changes; interp1(p_times, p_info.values, times)'];
    delays = [delays, p_info.delay];
end

% make the matrix X sucht that DeltaV = X * F
X = [];
for i=1:length(plasticity_files)
    R_d_T = zeros(size(gc_rates,2), length(learning_window_times) + 1);
    for j = 1:size(gc_rates,2)
        R_d_T(j,:) = [interp1(times - delays(i), gc_rates(:,j)', ...
            learning_window_times, 'nearest', 'extrap'), ...
            sum(gc_rates(:,j))]';
    end
    % TODO: Convolve with synaptic kernel
    X_i = gc_rates * R_d_T;
    X = cat(1, X, X_i);
end

pinvX = pinv(X);
learning_rule = pinvX * voltage_changes;
estimated_voltage_changes = X * learning_rule;

figure;
plot(learning_window_times, learning_rule(1:(length(learning_rule)-1)));