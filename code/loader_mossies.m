function [celltypes,rspstore]=loader_mossies(pth,avgtype)
folders=dir(pth);
folders = folders(find(cellfun(@isempty,strfind({folders.name},'.'))));
folders = folders(find(cellfun(@isempty,strfind({folders.name},'modeling'))));

%this is the time window over which we study the EOCD response:
winstart=-.025;
winend=.2;
dt=5e-5;

rspstore=cell(1);

cellcountnum=1;
for celltype=1:length(folders) %loop over folders

    files=dir([pth folders(celltype).name]);
    files(1:2)=[]; %remove current/parent directory from list
    files = files(find(cellfun(@isempty,strfind({files.name},'.DS_Store')))); %remove git artifacts
    
    
    for cellnum=1:length(files)
        
        %this code loads the data and figures out which channel is which.
        data=load([pth folders(celltype).name '/' files(cellnum).name]); %first load data as a struct to preserve channel names
        channelnames=fieldnames(data);
        for i=1:2
            channelnames{i}=str2double(regexprep(channelnames{i}(end-3:end),'\D','')); %trim nonnumerics to get the channel #
        end
        data=struct2cell(data); %now that we have param names, convert data to a cell so we can work with it
        if(channelnames{2}==10)%channel 10 contains spiking data, other channel is events
            eventind=1; spind=2;
        else
            eventind=2; spind=1;
        end

        eventtimes=data{eventind}.times;
        spiketimes=data{spind}.times;


        rspstore{cellcountnum}=zeros(length(eventtimes),ceil((winend-winstart)/dt));
        count=1;
        %the response to each recorded EOCD (within a window defined at the
        %top) is now stored as a row in a matrix; the response matrix for
        %each cell is stored within the cell array rspstore.
        for i=2:length(eventtimes)-1
            if((eventtimes(i)-eventtimes(i-1)>winend)&&(eventtimes(i+1)-eventtimes(i)>winend))  %only use well-isolated events
                mrkstart=find(spiketimes>eventtimes(i)+winstart,1);                             %find the first spike in the current EOCD window
                mrkend=find(spiketimes>eventtimes(i)+winend,1);                                 %find the first spike past the end of the current window
                spind=ceil((spiketimes(mrkstart:mrkend-1)-eventtimes(i)-winstart)/dt);          %find all spikes between those two
                rspstore{cellcountnum}(count,spind)=1;
                
                
                if(strcmp(avgtype,'alltrials')) %can decide in function call whether to store all trials or just trials where the MF spiked
                    count=count+1;
                elseif(~isempty(spind))
                    count=count+1;
                end
            end
        end
        ntrials=count-1;
        rspstore{cellcountnum}=rspstore{cellcountnum}(1:ntrials,:);
        
        %finally, format the name of this mossy fiber for display:
        menuname=folders(celltype).name;
        menuname(find(menuname=='_'):find(menuname=='_')+3)=[];
        celltypes{cellcountnum}=[menuname ' ' num2str(cellnum) ' ' files(cellnum).name(1:end-4)]; %append this to include MF numbers
        cellcountnum=cellcountnum+1;
    end
end