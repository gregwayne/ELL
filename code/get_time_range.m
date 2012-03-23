function times = get_time_range(filenames, method)

if strcmp(method, 'INTERSECTION')
    min_time = -Inf;
    max_time = Inf;
    interval = Inf;
    for i=1:length(filenames);
        data = importdata(char(filenames(i)));
        local_times = data.start + data.interval*((1:data.length)-1);
        
        min_time = max(min_time, local_times(1));
        max_time = min(max_time, local_times(length(local_times)));
        interval = min(interval, data.interval);
    end
    times = min_time:interval:max_time;
end

if strcmp(method, 'UNION')
    min_time = Inf;
    max_time = -Inf;
    interval = Inf;
    for i=1:length(filenames);
        data = importdata(char(filenames(i)));
        local_times = data.start + data.interval*((1:data.length)-1);
        
        min_time = min(min_time, local_times(1));
        max_time = max(max_time, local_times(length(local_times)));
        interval = min(interval, data.interval);
    end
    times = min_time:interval:max_time;
end
