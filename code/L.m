function loss = L(VS,PS,theta)

    parameters.v0       = theta(1);
    parameters.power    = theta(2);
    parameters.slope    = theta(3);
        
    loss = sum((PS-convert_voltages_to_spike_rates(VS,'POWER',parameters)).^2);

end