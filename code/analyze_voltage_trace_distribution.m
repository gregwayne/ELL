function analyze_voltage_trace_distribution(voltage_traces, times)

% Fourier transform
figure;
for i = 90:100 %10 %1:size(voltage_traces, 2)
    f = fft(voltage_traces(:,i)');
    n = ceil(length(f)/2);
    f = f(1:n);
    sampFreq = 1/(times(2) - times(1));
    freqs = sampFreq/2*linspace(0,1,n);
    plot(freqs(2:100), log(abs(f(2:100)).^2));
    hold all;
end
hold off;

%
for i = 1:size(voltage_traces, 2)
    voltage_traces(:,i) = voltage_traces(:,i) - mean(voltage_traces(1:100,i));
    %voltage_traces(:,i) = voltage_traces(:,i) / norm(voltage_traces(:,i));
end

[coef, score] = princomp(voltage_traces');

% plot the first few principle components
figure;
title('PCs');
for i=1:5
    plot(times, coef(:,i));
    hold all;
end
legend(num2str((1:5)'))
hold off;

% scatter plot in PC space
figure;
scatter(score(:,1)+score(:,2)+score(:,3), score(:,5),'filled');


% K-means
num_means = 6;
[IDX, C] = kmeans(voltage_traces', num_means, 'distance', 'cosine');
figure;
for i = 1:num_means
    plot(times, C(i,:),'DisplayName',num2str(i));
    hold all;
end
hold off;
legend(num2str((1:num_means)'))

figure;
hist(IDX);
