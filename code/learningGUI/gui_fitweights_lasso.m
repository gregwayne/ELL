function gui_fitweights_lasso(source,eventdata,uictrls,mean_mf,real_cells,mftypes)

disp(' ');
disp('fitting MFs/weights using Lasso...');

%get the latest model settings (and change some things for simulation)
GC_model=get_model_from_gui(uictrls);

DFmax=3;
[GC_model,~]=fitterlasso(GC_model,mean_mf,real_cells,DFmax);

%let the folks at home know what happened
disp(' ')
disp(['Fit weights for cell ' num2str(GC_model.GC_to_model) '!']);
disp('---------');
for i=1:3
    if(GC_model.MF_input(i)~=0)
        disp([mftypes{GC_model.MF_input(i)} ': ' num2str(GC_model.Ws(i))]);
    else
        disp('empty');
    end
end
disp('---------');

%update the sliders
set(0,'CurrentFigure',uictrls.fig);
for i=1:3
    set(uictrls.synapse_sliders{i},'Value',GC_model.Ws(i));
    set(uictrls.synapse_numbers{i},'Value',GC_model.MF_input(i)+1);
end