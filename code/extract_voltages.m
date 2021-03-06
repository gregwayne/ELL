function voltages = extract_voltages(filenames, times)
% Extract a matrix of voltages at each time, with one column per file

voltages = zeros(length(times), length(filenames));
for i=1:length(filenames)
    data = importdata(char(filenames(i)));
    tmp_times = data.start + data.interval*((1:data.length)-1);
    voltages(:,i) = interp1(tmp_times, data.values, times,'nearest','extrap');
    if strcmp(data.units, ' volt')
        voltages(:,i) = voltages(:,i) * 100; %these were actually in V/100
    elseif not(strcmp(data.units, 'mV'))
        data
    end
end
