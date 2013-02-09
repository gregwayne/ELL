function gc_rates = load_simulated_gcs()

gc_rates = [];
d = dir('../gcsWithTonic/*.mat');
for fIdx = 1:length(d)
    load(strcat('../gcsWithTonic/', d(fIdx).name));
    %for i = 1 : (size(meanSpikes,2) / 10)
    %    tmp_gc_rates(i,:) = sum(meanSpikes(:, (10 * (i-1) + 1) : (10 * i)),2);
    %end
    gc_rates = [gc_rates, meanSpikes'];
end
gc_rates = gc_rates(:, max(gc_rates) > 0); %get rid of cells thata don't spike
gc_rates = full(gc_rates);

end
