clear all;
close all;
raw_gc_filenames = {...
    '../rawGCdata/*.mat'};

% make list of files
raw_gc_files = get_file_list(raw_gc_filenames);

GCs = {};
SPs = {};

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
    
    GCs{i}      = raw_gc;
    
    % TESTING
    LT = 100000;
    figure(1);
    subplot(10,1,i);
    plot(raw_gc.lowgain.values(1:LT))
    hold on;
    
    spidxs = find_spike_indices(raw_gc.lowgain.values(1:LT));
    plot(spidxs,raw_gc.lowgain.values(spidxs),'.r');
    
    SPs{i}      = spidxs;
    
    % TODO
    % The EOD times are in raw_gc.cmdfilt.times
    % Need to find trial markers and minimum trial lengths
    % Then need to average traces over these trials
    % and fit probability of spiking based on the voltages
    
end

%% Maximum Likelihood
% TODO: Compute log-likelihood and gradient w.r.t. parameters
% Optimize convert_voltages_to_spike_rates.m

