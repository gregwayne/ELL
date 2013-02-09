clear all;
close all;
FULL = 1;
% granule_cell_filenames = {...
%     '../DATA2/CD_ipsp/*.mat',...
%     '../DATA2/CDpause_GCs/*.mat',...
%     '../DATA2/med/*.mat',...
%     '../DATA2/PCA/*.mat',...
%     '../DATA2/PCA_CDpause_GCs/*.mat',...
%     '../DATA2/PCA_late/*.mat',...
%     '../DATA2/PCA_med/*.mat',...
%     '../DATA2/PCA_med_late/*.mat'};
% %     'dataformodel/CDpause_GCs/*.mat',...
% %     'dataformodel/med_gcs/*.mat',...
% %     'dataformodel/PCA_CDpause_GCs/*.mat',...
% %     'dataformodel/PCA_GCs/*.mat',...
% %     'dataformodel/pca_med_gcs/*.mat',...
% %     'dataformodel/PCA_w_late/*.mat'};
% 
plasticity_filenames = {...
    '../dataformodel/plasticity/20111213/*.mat',...
    '../dataformodel/plasticity/20111220/*.mat'};

% % make lists of files
plasticity_files = get_file_list(plasticity_filenames);

% find the maximum and minimum times
times = get_time_range(plasticity_files, 'INTERSECTION');
% 
% % create matrices for granule cell voltages and firing rates
% gc_voltages = extract_voltages(gc_files, times);
% 
% % analyze_voltage_trace_distribution(gc_voltages, times);

gc_rates = load_simulated_gcs();

synaptic_kernel_parameters = struct('tau_on', 0.01, 'tau_off', 0.03);
%% fit each plasticity trace will least squares
%figure;
%title('Fits to each plasticity trace')
pinv_gc_rates = pinv(gc_rates);
delays = [];
voltage_changes = zeros(length(times),length(plasticity_files));
for i=1:length(plasticity_files)
    p_info = importdata(char(plasticity_files(i)));
    p_times = p_info.start + p_info.interval*((1:p_info.length)-1);
    
    voltage_changes(:,i) = interp1(p_times, p_info.values, times);
    delays = [delays, p_info.delay];

    %weight_changes = pinv_gc_rates * voltage_changes(:,i);
    %fitted_voltage_change =  gc_rates * weight_changes;
    %subplot(ceil(length(plasticity_files)/3), 3, i);
    %plot(times, voltage_changes(:,i), times, fitted_voltage_change);
    %title(p_info.delay);
end

unique_delays = sort(unique(delays));

%% fit the learning rule by least squares
learning_window_times = -0.10:0.001:0.10;
penalty = 2.0;

% construct vector by concatenating voltage changes at each delay

[learning_rule, predicted_voltage_changes] = ...
    fit_learning_rule_by_least_squares(gc_rates, times, ...
    voltage_changes, delays, learning_window_times, penalty, ...
    synaptic_kernel_parameters);

figure;
plot(learning_window_times, learning_rule(1:(length(learning_rule)-1))+...
    learning_rule(length(learning_rule)));
title('Plasticity rule with least squares')

%% plot resulting plasticity curves based on learned weights
figure;
title('Fits with least squares plasticity rule')
for i=1:length(plasticity_files)
    subplot(ceil(length(plasticity_files)/3), 3, i);
    plot(times, voltage_changes(:,i), ...
        times, predicted_voltage_changes(:,i));
    title(delays(i));
end

if FULL
%% test applying the plasticity rule
for whether_nonassociative_plasticity = 1 %0:1
    [learning_rule_parameters, predicted_voltage_changes] = fit_exp_learning_rule(...
        gc_rates, times, voltage_changes, delays, synaptic_kernel_parameters, ...
        whether_nonassociative_plasticity);
    figure;
    hold on;
    title('Optimal exponential learning rule predictions')
    for i=1:length(unique_delays)
        subplot(ceil(length(unique_delays)/2), 2, i);
        for delay_idx = find(delays == unique_delays(i))
            plot(times, voltage_changes(:,delay_idx),'b');, ...
            hold on;
            plot(times, predicted_voltage_changes(:,delay_idx),'g');
            hold on;
        end
        title(unique_delays(i));
    end
    subplot(ceil(length(unique_delays)/2), 2, i+1);
    learning_window_times = -0.1:(times(2)-times(1)):0.1;
    fit_exp_window = learning_rule_parameters.non_associative_weight + ...
        learning_rule_parameters.post_before_pre_amp * (learning_window_times - learning_rule_parameters.offset < 0) .* exp((learning_window_times - learning_rule_parameters.offset)/ learning_rule_parameters.tau_post_before_pre) + ...
        learning_rule_parameters.pre_before_post_amp * (learning_window_times - learning_rule_parameters.offset > 0) .* exp(-(learning_window_times - learning_rule_parameters.offset)/ learning_rule_parameters.tau_pre_before_post);
    plot(learning_window_times, fit_exp_window);
    hold on;
    plot(learning_window_times, 0);
    title('Learning window');
    xlabel('tpost - tpre');
    hold off;
    hgsave('full-fits.fig')
end
end
%% analyze stability of the learning rule with the given set of gc inputs
%stability_analysis(gc_rates, times, learning_rule_parameters);
