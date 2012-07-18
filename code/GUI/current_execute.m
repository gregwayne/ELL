GC_model_initialize();

if(~exist('mean_mf'))
    data_initialize();
end

make_gui;

v_fig = figure();
[realtrace,modeltrace,tran]=simulate_current_based_convolution(GC_model,mean_mf,real_cells);
set(0,'CurrentFigure',v_fig);
clf;
hold on;

h1=plot(tran,realtrace-mean(realtrace(1:200)));
h2=plot(tran,modeltrace-mean(modeltrace(1:200)),'g');

legend('real cell','fake cell','location','northeast');
xlabel('Time (s)');
ylabel('Vm (mV relative to resting)');
xlim([-.025 .15]);
% fixfig;
drawnow;

while 1
    if(quitsim)
       close all;
       break; 
    end
    for i=1:3
        GC_model.Ws(i)          = get(uictrls.synapse_sliders{i},'Value');
        GC_model.MF_input(i)    = get(uictrls.synapse_numbers{i},'Value')-1;
    end
    
    GC_model.tau_s = get(uictrls.tau_syn_slider,'Value');
%     GC_model.tau_s_ratio = get(uictrls.tau_syn_ratio,'Value');
    GC_model.tau_m = get(uictrls.tau_mem_slider,'Value');
    set(uictrls.tau_syn_label,'String',num2str(GC_model.tau_s));
    set(uictrls.tau_mem_label,'String',num2str(GC_model.tau_m));
    GC_model.spiking_on = get(uictrls.spikes,'Value');
    
    GC_model.GC_to_model=get(uictrls.GC_to_use,'Value');
    GC_model.nreps=str2num(get(uictrls.runs_to_avg,'String'));
    
    [realtrace,modeltrace,tran]=simulate_current_based_convolution(GC_model,mean_mf,real_cells);
    
    set(0,'CurrentFigure',v_fig);
    set(h1,'YData',realtrace-mean(realtrace(1:200)));
    set(h2,'YData',modeltrace-mean(modeltrace(1:200)));
    drawnow;
end
