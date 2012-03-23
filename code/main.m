clear all;
close all;
granule_cell_filenames = {...
    '../DATA2/CD_ipsp/*.mat',...
    '../DATA2/CDpause_GCs/*.mat',...
    '../DATA2/med/*.mat',...
    '../DATA2/PCA/*.mat',...
    '../DATA2/PCA_CDpause_GCs/*.mat',...
    '../DATA2/PCA_late/*.mat',...
    '../DATA2/PCA_med/*.mat',...
    '../DATA2/PCA_med_late/*.mat'};
%     'dataformodel/CDpause_GCs/*.mat',...
%     'dataformodel/med_gcs/*.mat',...
%     'dataformodel/PCA_CDpause_GCs/*.mat',...
%     'dataformodel/PCA_GCs/*.mat',...
%     'dataformodel/pca_med_gcs/*.mat',...
%     'dataformodel/PCA_w_late/*.mat'};

plasticity_filenames = {...
    '../dataformodel/plasticity/20111213/*.mat',...
    '../dataformodel/plasticity/20111220/*.mat'};

% make lists of files
gc_files = get_file_list(granule_cell_filenames);
plasticity_files = get_file_list(plasticity_filenames);

% find the maximum and minimum times
times = get_time_range(plasticity_files, 'INTERSECTION');

% create matrices for granule cell voltages and firing rates
gc_voltages = extract_voltages(gc_files, times);

%% convert to rates
gc_rates = zeros(size(gc_voltages));
%rate_parameters = struct('r0', 10, 'DeltaV', 5);
rate_parameters = struct('power',1,'v0',2, 'slope',1);

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
delays = [];
voltage_changes = zeros(length(times),length(plasticity_files));
for i=1:length(plasticity_files)
    p_info = importdata(char(plasticity_files(i)));
    p_times = p_info.start + p_info.interval*((1:p_info.length)-1);
    
    voltage_changes(:,i) = interp1(p_times, p_info.values, times);
    delays = [delays, p_info.delay];

    weight_changes = pinv_gc_rates * voltage_changes(:,i);
    fitted_voltage_change =  gc_rates * weight_changes;
    subplot(ceil(length(plasticity_files)/3), 3, i);
    plot(times, voltage_changes(:,i), times, fitted_voltage_change);
    title(p_info.delay);
end

%% fit the learning rule by least squares
learning_window_times = -0.05:0.005:0.05;
penalty = 10;

% construct vector by concatenating voltage changes at each delay

[learning_rule, predicted_voltage_changes] = ...
    fit_learning_rule_by_least_squares(gc_rates, times, ...
    voltage_changes, delays, learning_window_times, penalty);

figure;
plot(learning_window_times, learning_rule(1:(length(learning_rule)-1))+...
    learning_rule(length(learning_rule)));

%% plot resulting plasticity curves based on learned weights
figure;
for i=1:length(plasticity_files)
    
    subplot(ceil(length(plasticity_files)/3), 3, i);
    plot(times, voltage_changes(:,i), ...
        times, predicted_voltage_changes(:,i));

    title(delays(i));
end

%% test applying the plasticity rule

stdp_voltage_changes = zeros(size(voltage_changes));
weight_changes = zeros(size(gc_rates, 2), size(voltage_changes,2));
learning_rule_parameters = struct('potentiation_per_spike', 2e-4, 'amplitude', -1e-3, 'tau', 0.01)
figure;
title('Application of anti STDP rule')
for i=1:length(delays)
    weight_changes(:,i) = apply_plasticity_rule(times, gc_rates, delays(i), 'EXP_WINDOW', learning_rule_parameters);
    stdp_voltage_changes(:,i) = convolve_with_synaptic_kernel(gc_rates, times(2) - times(1), 0.003, 0.02) ...
        * weight_changes(:,i);
    subplot(ceil(length(plasticity_files)/3), 3, i);
    plot(times, voltage_changes(:,i), times, stdp_voltage_changes(:,i));
end
