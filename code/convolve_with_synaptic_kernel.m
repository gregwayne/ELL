function convolved_rates = convolve_with_synaptic_kernel(input_activity, dt, tau_fast, tau_slow)

% Construct synaptic kernel
filter_times = (-10*tau_slow):dt:(10*tau_slow);
synaptic_kernel = (filter_times > 0) .* (exp(-filter_times/tau_slow) - exp(-filter_times/tau_fast));
synaptic_kernel = synaptic_kernel / sum(synaptic_kernel);

% convolve firing rates with synaptic kernel
convolved_rates = zeros(size(input_activity));
for i=1:size(input_activity,2)
    convolved_rates(:,i) = conv(input_activity(:,i), synaptic_kernel, 'same');
end

end
