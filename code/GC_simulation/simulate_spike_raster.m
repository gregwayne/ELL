function raster = simulate_spike_raster(GC_model,rspstore,numtrials)

tmin = GC_model.min_t;
dt   = GC_model.dt;
tmax = GC_model.max_t;
tran = tmin+dt:dt:tmax;

raster=zeros(numtrials,length(tran));

for i=1:numtrials
    modeldat = simulate_current_based_expeuler(GC_model,rspstore);
    raster(i,find(modeldat==0))=1;
end