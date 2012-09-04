function raster = simulate_spike_raster(GC_model,rspstore,real_cells,numtrials)


raster=zeros(numtrials,size(real_cells,2));

for i=1:numtrials
    [~,modeldat,~] = simulate_current_based_expeuler(GC_model,rspstore,real_cells);
    raster(i,find(modeldat==0))=1;
end