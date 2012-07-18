GC_model={};
GC_model.tau_m       = 8.7; % 8.7 (ms)
GC_model.V_thresh    = -43; % (mV);
GC_model.rmgs        = 1; % unitless
GC_model.E_L         = -63; % (mV)
GC_model.V_reset     = GC_model.E_L; % (mV);
GC_model.tau_s       = 0.25; % (ms) -- look up
GC_model.Ws          = GC_model.rmgs*[.5 .5 .5];
GC_model.spiking_on  = 0;
GC_model.nreps       = 10;

GC_model.MF_input    = [0 0 0]; %0 for no input
GC_model.GC_to_model = 54;

if(~exist('rspstore'))
    data_initialize();
end

plotwin={};
plotwin.v_fig = figure();
[realtrace,modeltrace,tran]=simulate_conductance_based(GC_model,rspstore,real_cells);
set(0,'CurrentFigure',plotwin.v_fig);
clf;
hold on;

plotwin.h1=plot(tran,realtrace-mean(realtrace(1:200)));
plotwin.h2=plot(tran,modeltrace-mean(modeltrace(1:200)),'g');
legend('real cell','fake cell','location','northeast');
xlabel('Time (s)');
ylabel('Vm (mV relative to resting)');
xlim([-.025 .15]);
fixfig;
drawnow;

make_gui;




while 1
    if(quitsim)
%        clear all;
       close all;
       break; 
    end
    for i=1:3
        GC_model.Ws(i)          = get(uictrls.synapse_sliders{i},'Value');
        GC_model.MF_input(i)    = get(uictrls.synapse_numbers{i},'Value')-1;
    end
    
    GC_model.tau_s = get(uictrls.tau_syn_slider,'Value');
    GC_model.tau_m = get(uictrls.tau_mem_slider,'Value');
    set(uictrls.tau_syn_label,'String',num2str(GC_model.tau_s));
    set(uictrls.tau_mem_label,'String',num2str(GC_model.tau_m));
    GC_model.spiking_on = get(uictrls.spikes,'Value');
    
    GC_model.GC_to_model=get(uictrls.GC_to_use,'Value');
    GC_model.nreps=str2num(get(uictrls.runs_to_avg,'String'));
    
    [realtrace,modeltrace,tran]=simulate_conductance_based(GC_model,rspstore,real_cells);

    set(0,'CurrentFigure',plotwin.v_fig);
    set(plotwin.h1,'YData',realtrace-mean(realtrace(1:200)));
    set(plotwin.h2,'YData',modeltrace-mean(modeltrace(1:200)));
    drawnow;
end
