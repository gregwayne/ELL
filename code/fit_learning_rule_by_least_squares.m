function [learning_rule, predicted_voltage_changes] = ...
    fit_learning_rule_by_least_squares(gc_rates, gc_times, ...
    voltage_changes, delays, learning_window_times, penalty, ...
    synaptic_kernel_parameters)


num_plasticity_results = size(voltage_changes, 2);

% vecotrize the voltage changes
voltage_changes = reshape(voltage_changes, [], 1);

% Convolve rates with synaptic kernel
dt = gc_times(2) - gc_times(1);
convolved_rates = convolve_with_synaptic_kernel(gc_rates, dt, ...
    synaptic_kernel_parameters.tau_on, synaptic_kernel_parameters.tau_off);

% make the matrix X sucht that DeltaV = X * Learning_Rule
X = [];
for i=1:num_plasticity_results
    R_d_T = zeros(size(gc_rates,2), length(learning_window_times) + 1);
    for j = 1:size(gc_rates,2)
        R_d_T(j,:) = [interp1(gc_times - delays(i), ...
            gc_rates(:,j)', learning_window_times, 'nearest', 'extrap'), ...
            sum(gc_rates(:,j))]';
    end

    X_i = convolved_rates * R_d_T;
    X = cat(1, X, X_i);
end
breakTimes = linspace(-0.005, 0.005, 5);
err = Inf
for bIdx = 1:length(breakTimes)
    breakTime = breakTimes(bIdx);

    % penalize differences between adjacent terms
    for i=1:(length(learning_window_times)-1)
        if or(learning_window_times(i) > breakTime, learning_window_times(i+1) <= breakTime)
            X_i = zeros(1,length(learning_window_times) + 1);
            X_i(i) = penalty;
            X_i(i+1) = - penalty;
            X = cat(1, X, X_i);
        end
    end
    X_i = zeros(1,length(learning_window_times) + 1);
    X_i(end) = penalty;
    X_i(1) = - penalty;
    X = cat(1, X, X_i);

    %calculate learning rule
    pinvX = pinv(X);
    target = cat(1, voltage_changes, zeros(size(X,1) - size(voltage_changes,1), 1));
    tmpLearning_rule = pinvX * target;
        
    tmpErr = norm(X * tmpLearning_rule - target);
    if tmpErr < err
       learning_rule = tmpLearning_rule; 
       err = tmpErr;
    end

end

predicted_voltage_changes = ...
    reshape(X(1:length(voltage_changes),:) * learning_rule, [], ...
        num_plasticity_results);

end
