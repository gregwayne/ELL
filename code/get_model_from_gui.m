function GC_model=get_model_from_gui(uictrls)

GC_model={};
GC_model.V_thresh    = -43; % (mV);
GC_model.rmgs        = 1; % unitless
GC_model.E_L         = -63; % (mV)
GC_model.V_reset     = GC_model.E_L; % (mV);
for i=1:3
    GC_model.Ws(i)          = get(uictrls.synapse_sliders{i},'Value');
    GC_model.MF_input(i)    = get(uictrls.synapse_numbers{i},'Value')-1;
end
GC_model.tau_s = get(uictrls.tau_syn_slider,'Value');
GC_model.tau_m = get(uictrls.tau_mem_slider,'Value');
GC_model.GC_to_model=get(uictrls.GC_to_use,'Value');
GC_model.spiking_on = get(uictrls.spikes,'Value');
GC_model.nreps=str2num(get(uictrls.runs_to_avg,'String'));