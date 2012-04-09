function spike_rates = convert_voltages_to_spike_rates(...
    voltages, method, parameters)

if strcmp(method, 'THRESHOLD')
    vmin = min(voltages);
    vmax = max(voltages);
    thresh = vmin + parameters.threshold * (vmax - vmin);
    spike_rates = parameters.amplitude * (voltages > thresh);
end

if strcmp(method, 'EXP')
    spike_rates = parameters.r0 * exp(voltages / parameters.DeltaV);
end

if strcmp(method, 'POWER')
    spike_rates = (parameters.slope) * (voltages > parameters.v0) .* ...
        (voltages - parameters.v0) .^ parameters.power;
end

end
