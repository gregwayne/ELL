function gui_fitweights_grid(source,eventdata,uictrls,mean_mf,real_cells,mftypes)

disp(' ');
disp('fitting weights using grid search...');

%get the latest model settings (and change some things for simulation)
GC_model=get_model_from_gui(uictrls);

ncells=length(nonzeros(GC_model.MF_input));
wmin=get(uictrls.synapse_sliders{1},'Min')*ones(ncells,1);
wmax=get(uictrls.synapse_sliders{1},'max')*ones(ncells,1);
wbounds=[wmin(1) wmax(1)];

%future feature: restrict time window used for fitting?
% tstart=GC_model.min_t;
% tend=GC_model.max_t;
% disp(['using time window ' num2str(tstart/1e3) 's to ' num2str(tend/1e3) 's (relative to EOCD)']);
% tstart=(tstart-GC_model.min_t)/GC_model.dt+1;
% tend=(tend-GC_model.min_t)/GC_model.dt;

gridres=4;
nreps=8; %increase this to increase gridsearch resolution

%and now pass it on to gridfitter!
% GC_model=gridfitter(GC_model,mean_mf,real_cells,wmin,wmax,wbounds,gridres,nreps);
GC_model=MSEfitter(GC_model,mean_mf,real_cells);


%update the sliders
set(0,'CurrentFigure',uictrls.fig);
for i=1:3
    set(uictrls.synapse_sliders{i},'Value',GC_model.Ws(i));
end
disp('weights found! Final fit:');
disp('---------');
for i=1:3
    if(GC_model.MF_input(i)~=0)
        disp([mftypes{GC_model.MF_input(i)} ': ' num2str(GC_model.Ws(i))]);
    else
        disp('empty');
    end
end
disp(['With final grid resolution +/- ' num2str((wmax(1)-wmin(1))*(2/gridres)^nreps) '']);
disp('---------');