function stability_analysis(rates_in, times, learning_rule_parameters)

% discard rates that are always zero
rates = [];
for i = 1:size(rates_in,2)
    if sum(rates_in(:,i)) > 0
        rates = cat(2, rates, rates_in(:,i));
    end
end

num_neurons = size(rates,2);
dt = times(2) - times(1);

%% periodize the rates, epsp functions, and learning window
%zero_idx = find(times == 0);
%rates = rates(zero_idx:size(rates,1), :);
%times = times(zero_idx:length(times));

epsp_shape = (times > 0) .* (exp(-times / 0.02) - exp(-times / 0.002));
epsp_shape = epsp_shape / sum(epsp_shape) / dt;
plasticity_window = (times > 0) .* exp( - times / learning_rule_parameters.tau_pre_before_post) * learning_rule_parameters.pre_before_post_amp ...
    + (times < 0) .* exp(times / learning_rule_parameters.tau_post_before_pre) * learning_rule_parameters.post_before_pre_amp;

% Fourier transform rates,
ft_rates = fft(rates);
ft_epsp = fft(epsp_shape);
ft_plasticity = fft(plasticity_window);

freq = linspace(0, 1, ceil(length(ft_epsp/2))) /  dt;

% % take only low frequency elements
% MAX_FREQ = 100;
% MAX_IDX = min(find(freq > MAX_FREQ));
% freq = freq(1:MAX_IDX);
% ft_rates = ft_rates(1:MAX_IDX,:);
% ft_epsp = ft_epsp(1:MAX_IDX);
% ft_plasticity = ft_plasticity(1:MAX_IDX);

% construct matrix for the local derivative of the fourier coefficients of
% a perturbation of the voltage
deriv_matrix = diag(ft_epsp) * ft_rates * ctranspose(diag(ft_plasticity) * ft_rates);

[U,S,V] = svd(diag(ft_epsp) * ft_rates,'econ');
%[U,S,V] = svd(diag(ft_epsp) * ft_rates);
Q = ctranspose(U) * deriv_matrix * U;

[eigenvectors, eigenvalues] = eig(Q);

[rd, sort_idx] = sort(real(diag(eigenvalues)), 'ascend');

% sort the eigenvectors
sorted_eigenvectors = cat(1, eigenvectors, sort_idx');
sorted_eigenvectors = sortrows(sorted_eigenvectors', size(sorted_eigenvectors,1))';
sorted_eigenvectors = sorted_eigenvectors(1:(end-1),:);


figure; pcolor(abs(sorted_eigenvectors)); colorbar;

figure;
plot(rd, '.');

disp(max(rd));

figure;
plot(log(-real(rd))/log(10));
title('log_10(-real(eigenvalues))')

max_idx = min(find(freq > 1000));

figure;
power_spectrum = abs(U*sorted_eigenvectors).^2;
plot(freq(1:max_idx),log(power_spectrum(1:max_idx,end)),'r');
hold on;
plot(freq(1:max_idx),log(power_spectrum(1:max_idx,end-1)),'b');
plot(freq(1:max_idx),log(power_spectrum(1:max_idx,1)),'g');
plot(freq(1:max_idx),log(power_spectrum(1:max_idx,2)),'k');
hold off;