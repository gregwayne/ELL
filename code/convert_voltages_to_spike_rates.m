function spike_rates = convert_voltages_to_spike_rates(...
    voltages, method, parameters)

if strcmp(method, 'EXP')
    spike_rates = parameters.r0 * exp(voltages / parameters.DeltaV);
end

if strcmp(method, 'POWER')
    spike_rates = (voltages > parameters.v0) .* ...
        (voltages - parameters.v0) .^ parameters.power;
end

end