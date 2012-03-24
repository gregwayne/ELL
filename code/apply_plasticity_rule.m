function weight_changes = apply_plasticity_rule(times, rates, ...
    mg_spike_times, rule, parameters)

if strcmp(rule,'EXP_WINDOW')

    weight_changes = parameters.potentiation_per_spike * ...
        sum(rates,1) * (times(2) - times(1));
    
    for i = 1:length(mg_spike_times)
        t = mg_spike_times(i);
        weight_changes = weight_changes + parameters.amplitude * ...
            ((times < t) .* exp((times -t)/parameters.tau)) * rates * (times(2) - times(1));
    end
end
if strcmp(rule,'GENEREAL')

    weight_changes = parameters.potentiation_per_spike * ...
        sum(rates,1) * (times(2) - times(1));
    
    for i = 1:length(mg_spike_times)
        t = mg_spike_times(i);
        weight_changes = weight_changes + interp(parameters.window_times, parameters.window, times - t, 'linear', 0) * rates;
    end

end
