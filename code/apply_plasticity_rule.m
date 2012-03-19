function weight_changes = apply_plasticity_rule(times, rates, ...
    mg_spike_times, rule, parameters)

if strcmp(rule,'EXP_WINDOW')

    weight_changes = parameters.potention_per_spike * ...
        sum(rates,1) * (times(2) - times(1));
    
    for i = 1:length(mg_spike_times)
        t = mg_spike_times(i);
        weight_changes = weight changes + parameters.amplitude * ...
            ((times < t) .* exp((times -t)/parameters.tau)) * rates;
    end

end