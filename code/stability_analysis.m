function stability_analysis(rates, times_c)
close all;

%times_c = times;
%rates = [];
for i = 1:size(gc_rates,2)
    if sum(gc_rates(:,i)) > 0
        rates = cat(2, rates, gc_rates(:,i) / sum(gc_rates(:,i)) / dt);
    end
end

num_neurons = size(rates,2);
dt = times_c(2) - times_c(1);

% periodize the rates, epsp functions, and learning window
zero_idx = find(times_c == 0);
rates = rates(zero_idx:size(rates,1), :);
times_c = times_c(zero_idx:length(times_c));

epsp_shape = exp(-times_c / 0.02) - exp(-times_c / 0.002);
epsp_shape = epsp_shape / sum(epsp_shape) / dt;
plasticity_window = -epsp_shape;

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
Q = ctranspose(U) * deriv_matrix * U;

[eigenvectors, eigenvalues] = eig(Q);

[rd, sort_idx] = sort(real(diag(eigenvalues)), 'descend');

sorted_eigenvectors = sortrows(eigenvectors', sort_idx)';
figure; pcolor(abs(sorted_eigenvectors)); colorbar;

figure;
plot(rd, '.');

figure; plot(log(-real(eigenvalues))/log(10));
plot(log(-real(rd))/log(10));

