function [realtrace,modeltrace,tran] = simulate_current_based_expeuler(GC_model,rspstore,real_cells)

dt          = GC_model.dt;
min_t       = GC_model.min_t;
max_t       = GC_model.max_t;
times       = min_t+dt:dt:max_t;
nsteps      = length(times);
ninputs     = length(GC_model.MF_input);
inputs      = zeros(ninputs,nsteps);

tau_m       = GC_model.tau_m;
thr         = GC_model.V_thresh;

%make our inputs (convolved with synaptic kernel)
for i=1:ninputs
    if(GC_model.MF_input(i)>0)
        inputs(i,:) = draw_MF_input(rspstore{GC_model.MF_input(i)});
        inputs(i,:) = convolve_matrix_by_tau(dt,GC_model.tau_s,inputs(i,:));
    end
end

%simulate cell once, then check for spikes, and resimulate until we get
%them all
modeltrace = GC_model.E_L + GC_model.Ws*convolve_matrix_by_tau(dt,tau_m,inputs);

spiketimes=1;
count=2;
while(find(modeltrace>GC_model.V_thresh))
    
    tprevious=spiketimes(count-1);
    t=tprevious - 1 + find(modeltrace(tprevious:end) > thr,1);
    
    modeltrace(t)=-Inf;
    if(t+1<length(modeltrace))
        modeltrace(t+1:end) = GC_model.V_reset + GC_model.Ws*convolve_matrix_by_tau(dt,tau_m,inputs(:,t+1:end));
    end
    
    spiketimes(count)=t;
    count=count+1;
end

spiketimes(1)=[]; %get rid of the dummy spike
modeltrace(spiketimes)=0; %get rid of the -Infs

tran=times;
realtrace=real_cells(GC_model.GC_to_model,1:nsteps);
