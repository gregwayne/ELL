function plotPlasticityFits(source, eventdata, gc_rates, convolved_rates, ...
    voltage_changes, delays, times, plotHandle, learningRulePlot)

global uictrls;

learningRule = struct;
learningRule.tau_pre = str2num(get(uictrls.tau_pre_edit,'String'));
learningRule.tau_post = str2num(get(uictrls.tau_post_edit,'String'));
%learningRule.tau_pre_on = str2num(get(uictrls.tau_pre_on_edit,'String'));
%learningRule.tau_post_on = str2num(get(uictrls.tau_post_on_edit,'String'));
learningRule.amp_pre = str2num(get(uictrls.amp_pre_edit,'String'));
learningRule.amp_post = str2num(get(uictrls.amp_post_edit,'String'));
learningRule.offset = str2num(get(uictrls.offset_edit,'String'));
learningRule.amp_non_associative = str2num(get(uictrls.amp_non_edit,'String'));
    
offset = learningRule.offset;
tau_pre = learningRule.tau_pre;
tau_post = learningRule.tau_post;
tau_pre_on  = 0.005;
tau_post_on = 0.005;
amp_pre = learningRule.amp_pre;
amp_post = learningRule.amp_post;
amp_non_associative = learningRule.amp_non_associative;

unique_delays = sort(unique(delays));

longGCrates = [gc_rates(1:500,:); gc_rates];
longTimes = [times(1:500) + times(1) - times(501), times];
dt = times(2) - times(1);


get(uictrls.fitAmpSep, 'Value')
if get(uictrls.fitAmpSep, 'Value')
    set(0, 'CurrentFigure', learningRulePlot);
    clf;
    set(0, 'CurrentFigure', plotHandle);
    clf;
    for i=1:length(unique_delays)
        currentDelay = unique_delays(i);
        tmpVoltageChanges = voltage_changes(:, delays == currentDelay);

        num_plasticity_results = size(tmpVoltageChanges, 2)
        % make the matrix X sucht that DeltaV = X * Learning_Rule
        X = [];
        for traceIdx=1:num_plasticity_results
            post_minus_pre_delay = unique_delays(i) - longTimes;
            pre_kernel = (post_minus_pre_delay - offset < 0) .* (1 - exp((post_minus_pre_delay - offset) / tau_pre_on)) ...
                .* exp((post_minus_pre_delay - offset) / tau_pre)  * dt;
            post_kernel= (post_minus_pre_delay - offset > 0) .* (1 - exp((-post_minus_pre_delay+ offset) / tau_post_on)) ...
                .* exp((-post_minus_pre_delay+ offset) / tau_post) * dt;
            pre_kernel(isnan(pre_kernel)) = 0;
            post_kernel(isnan(post_kernel)) = 0;
            non_associative_kernel = ones(size(longTimes)) * dt;
            R_d_T = ([pre_kernel; post_kernel; non_associative_kernel] * longGCrates)';
            X_i = convolved_rates * R_d_T;
            X = cat(1, X, X_i); 
        end
        pinvX = pinv(X);
        learning_rule = pinvX * reshape(tmpVoltageChanges, [], 1);
        predicted_voltage_changes = ...
        reshape(X * learning_rule, [], num_plasticity_results);

        set(0, 'CurrentFigure', plotHandle);
        subplot(ceil(length(unique_delays)/2), 2, i);
        for idx = 1:num_plasticity_results
            plot(times, tmpVoltageChanges(:,idx),'b');, ...
            hold on;
            plot(times, predicted_voltage_changes(:,idx),'g');
            hold on;
        end

        set(0, 'CurrentFigure', learningRulePlot);
        subplot(ceil(length(unique_delays)/2), 2, i);
        learning_window_times = -0.1:(times(2)-times(1)):0.1;
        fit_exp_window = learning_rule(3) + ...
            + learning_rule(1)  * (learning_window_times - offset < 0) .* (1 - exp((learning_window_times - offset) / tau_pre_on)) ...
                .* exp((learning_window_times - offset) / tau_pre) ...
            + learning_rule(2) * (learning_window_times - offset > 0) .* (1 - exp((-learning_window_times+ offset) / tau_post_on)) ...
                .* exp((-learning_window_times+ offset) / tau_post);
        plot(learning_window_times, fit_exp_window);
        hold on;
        %plot(learning_window_times, 0);
        title('Learning window');
        xlabel('tpost - tpre');
        hold off;


        title(unique_delays(i));
    end

else
    num_plasticity_results = size(voltage_changes, 2);
    % make the matrix X sucht that DeltaV = X * Learning_Rule
    X = [];
    for i=1:num_plasticity_results

        post_minus_pre_delay = delays(i) - longTimes;

        pre_kernel = (post_minus_pre_delay - offset < 0) .* (1 - exp((post_minus_pre_delay - offset) / tau_pre_on)) ...
            .* exp((post_minus_pre_delay - offset) / tau_pre)  * dt;
        post_kernel= (post_minus_pre_delay - offset > 0) .* (1 - exp((-post_minus_pre_delay+ offset) / tau_post_on)) ...
            .* exp((-post_minus_pre_delay+ offset) / tau_post) * dt;
        pre_kernel(isnan(pre_kernel)) = 0;
        post_kernel(isnan(post_kernel)) = 0;
        non_associative_kernel = ones(size(longTimes)) * dt;
        R_d_T = ([pre_kernel; post_kernel; non_associative_kernel] * longGCrates)';
        X_i = convolved_rates * R_d_T;
        X = cat(1, X, X_i); 
    end

    if get(uictrls.fitAmpsGlobal,'Value')
        pinvX = pinv(X);
        learning_rule = pinvX * reshape(voltage_changes, [], 1);
        predicted_voltage_changes = ...
        reshape(X * learning_rule, [], num_plasticity_results);

        set(uictrls.amp_pre_edit, 'String', num2str(learning_rule(1)));
        set(uictrls.amp_post_edit, 'String', num2str(learning_rule(2)));
        set(uictrls.amp_non_edit, 'String', num2str(learning_rule(3)));
     
    else
        predicted_voltage_changes = reshape(X * [amp_pre; amp_post; amp_non_associative], [], num_plasticity_results);
    end

    set(0, 'CurrentFigure', plotHandle);
    clf;
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

end


end
