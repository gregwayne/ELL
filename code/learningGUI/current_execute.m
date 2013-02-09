addpath('..');
gc_rates = load_simulated_gcs();

plasticity_filenames = {...
    '../dataformodel/plasticity/20111213/*.mat',...
    '../dataformodel/plasticity/20111220/*.mat'};

% % make lists of files
plasticity_files = get_file_list(plasticity_filenames);

times = get_time_range(plasticity_files, 'INTERSECTION');
dt = times(2) - times(1);

delays = [];
voltage_changes = zeros(length(times),length(plasticity_files));
for i=1:length(plasticity_files)
    p_info = importdata(char(plasticity_files(i)));
    p_times = p_info.start + p_info.interval*((1:p_info.length)-1);
    
    voltage_changes(:,i) = interp1(p_times, p_info.values, times);
    delays = [delays, p_info.delay];
end

synaptic_kernel_parameters = struct('tau_on', 0.005, 'tau_off', 0.01);

convolved_rates = convolve_with_synaptic_kernel(gc_rates, dt, ...
    synaptic_kernel_parameters.tau_on, synaptic_kernel_parameters.tau_off);


plotHandle = figure();
learningRulePlot = figure();
make_gui;

%[realtrace,modeltrace,tran]=simulate_current_based_convolution(GC_model,mean_mf,real_cells);
%set(0,'CurrentFigure',v_fig);
%clf;
%hold on;

%h1=plot(tran,realtrace-mean(realtrace(1:200)));
%h2=plot(tran,modeltrace-mean(modeltrace(1:200)),'g');

legend('real cell','fake cell','location','northeast');
xlabel('Time (s)');
ylabel('Vm (mV relative to resting)');
xlim([-.025 .15]);
% fixfig;
% drawnow;

while 1
    if(quitsim)
       close all;
       break; 
    end
    %set(0,'CurrentFigure',v_fig);
    %set(h1,'YData',realtrace-mean(realtrace(1:200)));
    %set(h2,'YData',modeltrace-mean(modeltrace(1:200)));
    drawnow;
end
