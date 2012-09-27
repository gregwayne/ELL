function rspstore = replace_20120419_003(avgtype) %this is stupid

pth = '/Users/Ann/Downloads/20120419_003/';
subfiles = cleandir(pth);
subfiles = {subfiles.name};

winstart=-.025;
winend=.2;
dt=5e-5;


rspstore=zeros(1,4500);
count=1;
for input=1:length(subfiles)
    data=load([pth subfiles{input}]); %first load data as a struct to preserve channel names
    channelnames=fieldnames(data);

    for i=1:3
        channels(i)=str2double(regexprep(channelnames{i}(end-3:end),'\D','')); %trim nonnumerics to get the channel #
    end
    eventind = find(channels == 10);
    spind = find(channels==60);


    eventtimes = data.(channelnames{eventind}).times;
    spiketimes = data.(channelnames{spind}).times;


    for i=2:length(eventtimes)-1
        if((eventtimes(i)-eventtimes(i-1)>winend)&&(eventtimes(i+1)-eventtimes(i)>winend))  %only use well-isolated events
            mrkstart=find(spiketimes>eventtimes(i)+winstart,1);                             %find the first spike in the current EOCD window
            mrkend=find(spiketimes>eventtimes(i)+winend,1);                                 %find the first spike past the end of the current window
            spind=ceil((spiketimes(mrkstart:mrkend-1)-eventtimes(i)-winstart)/dt);          %find all spikes between those two
            rspstore(count,spind)=1;


            if(strcmp(avgtype,'alltrials')) %can decide in function call whether to store all trials or just trials where the MF spiked
                count=count+1;
            elseif(~isempty(spind))
                count=count+1;
            end
        end
    end

end