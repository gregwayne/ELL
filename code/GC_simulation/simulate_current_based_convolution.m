function [realtrace,modeltrace,tran] = simulate_current_based_convolution(GC_model,mean_mf,real_cells)

dt          = GC_model.dt;
min_t       = GC_model.min_t;
max_t       = GC_model.max_t;
times       = min_t:dt:max_t;
nsteps      = length(times)-1;
ninputs     = length(GC_model.MF_input);
inputs      = zeros(ninputs,nsteps);

vmean = zeros(1,nsteps);

realtrace = real_cells(GC_model.GC_to_model,1:nsteps);

for i=1:ninputs
    if(GC_model.MF_input(i)>0)
        inputs(i,:) = mean_mf(GC_model.MF_input(i),1:nsteps);
        inputs(i,:) = convolve_mossies(GC_model,inputs(i,:));
    end
end

vmean = GC_model.E_L + GC_model.Ws*inputs;

modeltrace = vmean(1:nsteps);
tran = (min_t+dt:dt:max_t)*10^-3;
