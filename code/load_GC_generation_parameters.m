function GCparams = load_GC_generation_parameters()

load '../GC_fitting_output/sept10_refit.mat';

theta = fit_mixing_coefficients_4_cell_types(Wsparse,mftypes);
MFindices = get_struct_of_celltypes(mftypes);
cellcats = fieldnames(MFindices);

% fit a gamma distribution to the fit weights
for i=1:length(cellcats)
    phat = gamfit(nonzeros(Wsparse(:,MFindices.(cellcats{i}))));
    gamma_k.(cellcats{i}) = phat(1);
    gamma_th.(cellcats{i}) = phat(2);
end

% and get the threshold statistics (this is slow and memory-intensive)
% [thr_mean,thr_var] = get_threshold_from_raw_GCs('../raw_gcs_mat');

thr_mean = 63 - 43.5;
thr_var = 6.8^2 + 4.9^2; %numbers from nate, stdev of baseline and threshold respectively


%these are all the parameters we need to make a new GC! 
GCparams.W_gamma_k    = gamma_k;
GCparams.W_gamma_th   = gamma_th;
GCparams.theta        = theta;
GCparams.thresh_mean  = thr_mean;
GCparams.thresh_var   = thr_var;