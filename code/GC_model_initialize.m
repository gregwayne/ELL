%this gets called when we launch the gui or the batch fitter, so we can
%load all the data and initialize the model. if you want to make changes to
%the default model settings, make them here!

gui_initialize();

GC_model={};
GC_model.tau_m       = 8.7; % 8.7 (ms)
GC_model.V_thresh    = -43; % (mV);
GC_model.rmgs        = 1; % unitless
GC_model.E_L         = -63; % (mV)
GC_model.V_reset     = GC_model.E_L; % (mV);
GC_model.tau_s       = 1; % (ms) -- look up
GC_model.Ws          = GC_model.rmgs*[.5 .5 .5];
GC_model.spiking_on  = 0;

GC_model.MF_input    = [0 0 0];
GC_model.GC_to_model = 1;

GC_model.dt          = 5e-2;
GC_model.min_t       = -.025e3;
GC_model.max_t       = .2e3;