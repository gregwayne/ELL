% load all the relevant data

data_initialize;                        %loads GC/MF data from files (this is really slow, boo)

GC_model_initialize;                    %sets up the GC_model struct we use for simulation

load ../GC_fitting_output/aug15.mat;    %load some saved fits of MF->GC weights.
                                        %aug15 is a pretty good recent
                                        %dataset; will update with the
                                        %final fits when they come.

time = size(real_cells,2);

[thr_mean,thr_var] = get_threshold_from_raw_GCs('../raw_gcs_mat'); %finally, get some info about our distribution of thresholds from the raw GC data


%% and generate bases! 

repeats_per_GC      = 5;                % we simulate the GC several times with random thresholds to increase our sample size.
trials_to_average   = 100;              % and we average spiking responses over many trials to get basis functions

bfns=zeros(numGCs*repeats_per_GC,time);

for cellnum=1:numGCs
    GC_model = load_weights_from_matrix(GC_model,Wsparse,cellnum);

    for rep=1:repeats_per_GC
        GC_model.V_thresh = GC_model.E_L + max(thr_mean + randn(1)*sqrt(thr_var),5); %don't let it pick a relative threshold of less than 5
        
        raster = simulate_spike_raster(GC_model,rspstore,real_cells,trials_to_average);
        bfns((cellnum-1)*repeats_per_GC+rep,:) = mean(raster);
    end
end