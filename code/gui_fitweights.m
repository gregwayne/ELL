function GC_model=GC_fitweights(source,eventdata,uictrls,rspstore,real_cells,plotwin)

%get the latest model settings (and change some things for simulation)
GC_model=get_model_from_gui(uictrls);
GC_model.spiking_on=0;
%I drop the number of reps you average over down to 5 for fitting, to save
%time.  if cells are particularly noisy this may be a bad idea.
GC_model.nreps=5;

ncells=length(nonzeros(GC_model.MF_input));
wmin=get(uictrls.synapse_sliders{1},'Min')*ones(ncells,1);
wmax=get(uictrls.synapse_sliders{1},'max')*ones(ncells,1);

%fit the weights via a nested grid search
% draws a coarse grid, finds the optimum point, then zooms in around that
% point, draws a new grid, and repeats.
% the last two arguments of optWeights are gridres and nreps- the number of
% points along each dimension we draw in our grid, and the number of
% iterations we perform.
% computing time should scale like nreps*(gridres^ncells), where ncells is
% the number of mossy fiber inputs (so, 1-3). I think the resolution of the
% final selected weights should be +/- (wmax-wmin)(2/gridres)^nreps-- so eg
% for wmax=10, wmin=0, gridres=4, and nreps=9, we get weights +/- 0.02.
% Probably there's a cute optimization problem in picking gridres/nreps to
% get the best reps/error tradeoff for a given number of inputs.
tic
GC_model=optWeights(uictrls,GC_model,rspstore,real_cells,plotwin,wmin,wmax,4,5);
toc
set(0,'CurrentFigure',plotwin.v_fig);
title('');

%update the sliders
set(0,'CurrentFigure',uictrls.fig);
for i=1:3
    set(uictrls.synapse_sliders{i},'Value',GC_model.Ws(i));
end