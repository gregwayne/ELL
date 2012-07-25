function spike_train = draw_mossy(mossy_pars)

    %          num_pulses: [1 double]
    %                  t0: [1 double]
    %                  e0: [1 double]
    %               delta: [1 double]
    %        jitter_stdev: [1 double]
    %             e_slope: [1 double]

    dt              = 5e-5; % (s)
    spike_train     = zeros(4500,1);
    np              = mossy_pars.num_pulses;
    t0              = 0.025 + mossy_pars.t0; % when EOCD occurs   
    delta           = mossy_pars.delta;
    e0              = mossy_pars.e0;
    e_slope         = mossy_pars.e_slope;
    jitter_stdev    = mossy_pars.jitter_stdev;
        
    for n=1:np
       
        % ask Patrick about this
        if e0 < 1
            ns = rand < e0;
        else
            ns  = poissrnd(e0);
        end
        
        for k=1:ns
           
            jt                          = jitter_stdev*randn;
            ts                          = t0 + jt;
            idx                         = max(1,ceil(ts/dt));
            spike_train(idx)            = 1;
                        
        end
        
        e0                          = e0 - e_slope;
        t0                          = t0 + delta;
        
    end
    
end