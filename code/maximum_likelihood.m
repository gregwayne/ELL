clear all;
close all;
raw_gc_filenames = {...
    '../rawGCdata/*.mat'};

% make list of files
raw_gc_files = get_file_list(raw_gc_filenames);

GCs     = {};
SPs     = {};
EODs    = {}; 

for i=1:length(raw_gc_files)
   
    raw_gc_info     = importdata(char(raw_gc_files(i)));
    raw_gc_fields   = fieldnames(raw_gc_info);

    raw_gc          = struct;
    a               = getfield(raw_gc_info,char(raw_gc_fields(1)));
    b               = getfield(raw_gc_info,char(raw_gc_fields(2)));
    if strcmp(a.title,'lowgain')
        raw_gc.lowgain = a;
        raw_gc.cmdfilt = b;
    else
        raw_gc.lowgain = b;
        raw_gc.cmdfilt = a;
    end
        
    % Fix format of cells recorded before 2010
    if strncmp(char(raw_gc_fields(1)),'V200',4)  % i.e., before 2010
        % Nate says the scale was actually mV/100
        raw_gc.lowgain.values   = raw_gc.lowgain.values*10^2; 
        raw_gc.lowgain.units    = 'mV';
    end
    
    GCs{i}              = raw_gc;     
    [spidxs,spikes]     = find_spike_indices(raw_gc.lowgain.values);    
    SPs{i}              = spidxs;
    
    % Find trial markers
    EOD_diffs   = raw_gc.cmdfilt.times(2:end) ...
                            - raw_gc.cmdfilt.times(1:(end-1));
    EOD_diffs(end+1) = raw_gc.lowgain.start ...
                           + raw_gc.lowgain.length*raw_gc.lowgain.interval; 
    
    min_diff    = min(EOD_diffs);
    idx_diff    = floor(min_diff/raw_gc.lowgain.interval)-1;
    if idx_diff < 1
        error('Trial with no separation.');
    end
    
    EOD_idxs    = floor((raw_gc.cmdfilt.times/raw_gc.lowgain.interval) ...
                                                -raw_gc.lowgain.start)+1;
    EODs{i}     = EOD_idxs;
    
    av_voltage_trace    = zeros(idx_diff,1);
    summed_spikes       = zeros(idx_diff,1);
    NT                  = length(raw_gc.cmdfilt.times);
    for j=1:(NT-1) % bothered by edge effect -- so cutting off last sample
        id                  = EOD_idxs(j);
        av_voltage_trace    = av_voltage_trace ...
                               + raw_gc.lowgain.values(id:(id+idx_diff-1));
        summed_spikes       = summed_spikes + spikes(id:(id+idx_diff-1));
    end
        
    av_voltage_trace        = av_voltage_trace/(NT-1);
    summed_spikes           = summed_spikes/(NT-1);

    % Note: Ann says that some cells have unrecorded EODs
    
    %% Use Bayes' Rule to get p(spike|voltage) histogram
    % p(spike|voltage) = p(voltage|spike)*p(spike)/p(voltage)
    subplot(10,1,i);
    pos_idxs                = find(summed_spikes>0);
    [N,X]                   = hist(av_voltage_trace,10);
    p_spike                 = sum(summed_spikes)/length(summed_spikes);
    W                       = N/sum(N); % normalized number in bin
    v_given_spikes          = av_voltage_trace(pos_idxs);
    [M]                     = hist(v_given_spikes,X);
    Q                       = M/sum(M);
    plot(X,p_spike*Q./W);
    axis([-80,-20,0,2e-2]);
    
end

%% Maximum Likelihood
% TODO: Compute log-likelihood and gradient w.r.t. parameters
% Optimize convert_voltages_to_spike_rates.m
