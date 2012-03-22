function [learning_rule, predicted_voltage_changes] = ...
    fit_learning_rule_by_least_squares(gc_rates, gc_times, ...
    voltage_changes, delays, learning_window_times, penalty)


num_plasticity_results = size(voltage_changes, 2);

% vecotrize the voltage changes
voltage_changes = reshape(voltage_changes, [], 1);

% Convolve rates with synaptic kernel
tau_fast = 0.003;
tau_slow = 0.02;
dt = gc_times(2) - gc_times(1);
convolved_rates = convolve_with_synaptic_kernel(gc_rates, dt, tau_fast, tau_slow);

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

% penalize differences between adjacent terms
for i=1:(length(learning_window_times)-1)
    X_i = zeros(1,length(learning_window_times) + 1);
    X_i(i) = penalty;
    X_i(i+1) = - penalty;
    X = cat(1, X, X_i);
end

pinvX = pinv(X);
learning_rule = pinvX * ...
    cat(1, voltage_changes, zeros(length(learning_window_times)-1, 1));
predicted_voltage_changes = ...
    reshape(X(1:length(voltage_changes),:) * learning_rule, [], ...
        num_plasticity_results);

end
