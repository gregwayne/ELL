function [GC_model,stats]=fitterlasso(GC_model,mean_mf,real_cells,DFmax)

H=convolve_mossies(GC_model,mean_mf);
H=bsxfun(@minus,H,mean(H,2));
dat=real_cells(GC_model.GC_to_model,:);

%lasso does the heavy lifting
[w,stats]=mylasso(H',dat,'DFmax',DFmax,'NumLambda',150,'LambdaRatio',5e-5);
stats.w=w; %save these too!

nonzerocells=find(w(:,1))';
newWs=w(nonzerocells,1)';
nonzerocells=[nonzerocells zeros(1,3-length(nonzerocells))];
newWs=[newWs zeros(1,3-length(newWs))];

GC_model.MF_input=nonzerocells;
GC_model.Ws=newWs;