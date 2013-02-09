gc_rates = load_simulated_gcs();
smooth_gc_rates = smoothts(gc_rates', 'b', 30)';

plasticity_filenames = {...
    '../dataformodel/plasticity/20111213/*.mat',...
    '../dataformodel/plasticity/20111220/*.mat'};
plasticity_files = get_file_list(plasticity_filenames);
times = get_time_range(plasticity_files, 'INTERSECTION');

figure;
subplot(221)
imagesc(times, times, log(smooth_gc_rates * smooth_gc_rates'));
axis equal
colorbar()

dt = times(2) - times(1)
synaptic_kernel_parameters = struct('tau_on', 0.01, 'tau_off', 0.03);
convolved_rates = convolve_with_synaptic_kernel(gc_rates, dt, ...
    synaptic_kernel_parameters.tau_on, synaptic_kernel_parameters.tau_off);

subplot(222)
imagesc(times, times, log(smooth_gc_rates * convolved_rates'));
axis equal
colorbar()

subplot(223)
imagesc(times, times, log(convolved_rates * convolved_rates'));
axis equal
colorbar()

subplot(224)
plot(times, mean(smooth_gc_rates,2))
hold all;
plot(times, mean(convolved_rates,2))
hold off;
