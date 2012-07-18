function gui_fitweights_LS(source,eventdata,uictrls,mean_mf,real_cells,mftypes)

disp(' ');
disp('fitting weights using least squares...');

%get the latest model settings (and change some things for simulation)
GC_model=get_model_from_gui(uictrls);


%and now pass it on to LS fitter!
GC_model=fitterLS(GC_model,mean_mf,real_cells);

MSE=compute_model_error(GC_model,mean_mf,real_cells,'MSE');

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
disp('---------');
disp(['MSE = ' num2str(MSE)]);