%load your GC generation parameters and your mossy fiber data:
load ../GC_fitting_output/GC_sim_params;
[mftypes, raw_mossy_traces, mean_mossy_traces, MF_indices] = mf_initialize();

%make a new fake GC, and then use it to simulate a spike raster
num_trials_to_simulate  = 50;
GC_model                = make_fake_GC(GCparams, MF_indices);
raster                  = simulate_spike_raster(GC_model, raw_mossy_traces, num_trials_to_simulate); 

%% If you'd like to look at the cell you generated:

%this code gives the average GC response
[~,GC_mean_trace,~] = simulate_current_based_convolution(GC_model,mean_mossy_traces,[]);

%and here's some figures. Hopefully raster occasionally has spikes in it!
figure(1);clf;
subplot(2,1,1);
plot(GC_mean_trace);
subplot(2,1,2);
imagesc(smoothts(raster,'e',100));






% You can poke around in the GC_model struct to get a sense of the cell you
% generate-- if you want to get a cell that spikes, try setting:
%
% GC_model.Ws       = [200 0 0]; %about a 20mV epsp
% GC_model.MF_input = [50 0 0];  %input from one of the early MF's, with the
%                                 other two channels empty (0 = no input)
% GC_model.V_thresh = -43;       %a reasonable threshold


% The GCparams stuct has the parameters we use to generate new GCs-- if you
% want to make the average threshold lower, change GCparams.thresh_mean
% (and probably also GCparams.thresh_var.)


% The code to generate spike rasters gets slower when the cell spikes a lot
% (because of the stupid way I chose to simulate it)-- because the cells
% are so sparse this hasn't been a problem for me, but if you have trouble
% with it let me know and I'll send you a better version.