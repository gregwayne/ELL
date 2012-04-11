function [learning_rule_parameters, predicted_voltage_changes] = ...
    fit_exp_learning_rule(gc_rates, gc_times, voltage_changes, delays, synaptic_kernel_parameters,...
    whether_nonassociative_plasticity)

% TODO: Option to get rid of non-associative component
% TODO: Option to shift discontinuity from zero


num_plasticity_results = size(voltage_changes, 2);

% vecotrize the voltage changes
voltage_changes = reshape(voltage_changes, [], 1);

% Convolve rates with synaptic kernel
dt = gc_times(2) - gc_times(1);
convolved_rates = convolve_with_synaptic_kernel(gc_rates, dt, ...
    synaptic_kernel_parameters.tau_fast, synaptic_kernel_parameters.tau_slow);


% Minimize over tau_pre and tau_post
objective_function = @(tau) ls_err(tau, gc_rates, gc_times, ...
    voltage_changes, delays, convolved_rates, whether_nonassociative_plasticity);

% Global grid search
tau_pre_values = 10.^(linspace(-3,1,15));
tau_post_values = 10.^(linspace(-3,1,15));
min_achieved_err = Inf;
for i=1:length(tau_pre_values)
    for j=1:length(tau_post_values)
        tmp_tau = [tau_pre_values(i), tau_post_values(j)];
        tmp_err = objective_function(tmp_tau);
        if tmp_err < min_achieved_err
            min_achieved_err = tmp_err;
            tau = tmp_tau;
        end
    end
end

% Local minimization
tau = fminunc(objective_function, tau); %simulannealbnd(objective_function, [2e-2, 2e-2]);

% calculate predicted voltages for optimal tau_pre and tau_post
tau_pre = tau(1);
tau_post = tau(2);

% make the matrix X sucht that DeltaV = X * Learning_Rule
X = [];
for i=1:num_plasticity_results

    post_minus_pre_delay = delays(i) - gc_times;

    pre_kernel = (post_minus_pre_delay < 0) .* exp((post_minus_pre_delay) / tau_pre)  * dt;
    post_kernel= (post_minus_pre_delay > 0) .* exp((-post_minus_pre_delay) / tau_post) * dt;

    if whether_nonassociative_plasticity
        non_associative_kernel = ones(size(gc_times)) * dt;
        R_d_T = ([pre_kernel; post_kernel; non_associative_kernel] * gc_rates)';
    else
        R_d_T = ([pre_kernel; post_kernel] * gc_rates)';
    end


    X_i = convolved_rates * R_d_T;
    X = cat(1, X, X_i); 
end

pinvX = pinv(X);
learning_rule = pinvX * voltage_changes;
predicted_voltage_changes = ...
    reshape(X * learning_rule, [], num_plasticity_results);
    

if whether_nonassociative_plasticity
    learning_rule_parameters = struct('tau_post_before_pre', tau_pre, 'tau_pre_before_post', tau_post, ...
        'post_before_pre_amp', learning_rule(1), 'pre_before_post_amp', learning_rule(2), ...
        'non_associative_weight', learning_rule(3));
else 
    learning_rule_parameters = struct('tau_post_before_pre', tau_pre, 'tau_pre_before_post', tau_post, ...
        'post_before_pre_amp', learning_rule(1), 'pre_before_post_amp', learning_rule(2), ...
        'non_associative_weight', 0);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OBJECTIVE FUNCTION TO BE MINIMIZED %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sq_err = ls_err(tau, gc_rates, gc_times, voltage_changes, ...
    delays, convolved_rates, whether_nonassociative_plasticity)

tau_pre = tau(1);
tau_post = tau(2);

dt = gc_times(2) - gc_times(1);

% make the matrix X sucht that DeltaV = X * Learning_Rule
X = [];
for i=1:length(delays)

    post_minus_pre_delay = delays(i) - gc_times;

    pre_kernel = (post_minus_pre_delay < 0) .* exp((post_minus_pre_delay) / tau_pre)  * dt;
    post_kernel= (post_minus_pre_delay > 0) .* exp((-post_minus_pre_delay) / tau_post) * dt;

    if whether_nonassociative_plasticity
        non_associative_kernel = ones(size(gc_times)) * dt;
        R_d_T = ([pre_kernel; post_kernel; non_associative_kernel] * gc_rates)';
    else
        R_d_T = ([pre_kernel; post_kernel] * gc_rates)';
    end

    X_i = convolved_rates * R_d_T;
    X = cat(1, X, X_i); 
end

pinvX = pinv(X);
learning_rule = pinvX * voltage_changes;
predicted_voltage_changes = X * learning_rule;
    

err = predicted_voltage_changes - voltage_changes;

sq_err = err' * err;

end 
