function GC_model=lassofitter(GC_model,mean_mf,real_cells,DFmax)

H=convolve_mossies(GC_model,mean_mf);
H=bsxfun(@minus,H,mean(H,2));
dat=real_cells(GC_model.GC_to_model,:);

%lasso does the heavy lifting
[w,stats]=lasso(H',dat,'DFmax',DFmax);

nonzerocells=find(w(:,1))';
newWs=w(nonzerocells,1)';
nonzerocells=[nonzerocells zeros(1,3-length(nonzerocells))];
newWs=[newWs zeros(1,3-length(newWs))];

GC_model.MF_input=nonzerocells;
GC_model.Ws=newWs;