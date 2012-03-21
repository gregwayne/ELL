function idxs = find_spike_indices(voltage_trace)

    LV                  = length(voltage_trace);
    
    derivative_filter   = [-1,1];
    derivatives         = conv(voltage_trace,derivative_filter);
    % find zero crossings for extrema
    ex_mask             = zeros(LV,1);
    for i=2:LV
        ex_mask(i)      = derivatives(i-1) < 0 & derivatives(i) >= 0;
    end
    % find ones that are local maxima
    second_derivatives  = conv(derivatives,derivative_filter);
    neg_mask            = second_derivatives < 0;
    max_mask            = ex_mask(1:LV) & neg_mask(1:LV);
    high_mask           = voltage_trace > -20; % spike threshold guess
        
    spikes              = max_mask & high_mask;
    idxs                = find(spikes > 0);
            
end