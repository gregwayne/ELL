function [celltypes,rspstore]=export_mossy_rasters(pth)
folders=dir(pth);
folders = folders(find(cellfun(@isempty,strfind({folders.name},'.'))));
folders = folders(find(cellfun(@isempty,strfind({folders.name},'modeling'))));

winstart=-.025;
winend=.2;
dt=5e-5;

rspstore=cell(94,1); %hacks!

cellcountnum=1;
for celltype=1:length(folders)
    files=dir([pth folders(celltype).name]);
    files(1:2)=[];
    files = files(find(cellfun(@isempty,strfind({files.name},'.DS_Store')))); %what why agh
    for cellnum=1:length(files)
        data=load([pth folders(celltype).name '/' files(cellnum).name]);
        channelnames=fieldnames(data);
        for i=1:2
            channelnames{i}=str2double(regexprep(channelnames{i}(end-3:end),'\D',''));%trim nonnumerics
        end
        data=struct2cell(data);


        if(channelnames{2}==10)%channel 10 contains spiking data, other channel is events
            eventind=1; spind=2;
        else
            eventind=2; spind=1;
        end

        eventtimes=data{eventind}.times;
        spiketimes=data{spind}.times;


        meanrsp=zeros(ceil((winend-winstart)/dt),1);
        rspstore{cellcountnum}=zeros(length(eventtimes),ceil((winend-winstart)/dt));
        count=1;
        for i=2:length(eventtimes)-1
            if((eventtimes(i)-eventtimes(i-1)>winend)&&(eventtimes(i+1)-eventtimes(i)>winend)) %only use well-isolated events
                mrkstart=find(spiketimes>eventtimes(i)+winstart,1);
                mrkend=find(spiketimes>eventtimes(i)+winend,1); %find the first spike past the end of the interval
                spind=ceil((spiketimes(mrkstart:mrkend-1)-eventtimes(i)-winstart)/dt);
                rspstore{cellcountnum}(count,spind)=1;
                count=count+1;
            end
        end
        ntrials=count-1;
        rspstore{cellcountnum}=rspstore{cellcountnum}(1:ntrials,:);
        menuname=folders(celltype).name;
        menuname(find(menuname=='_'):find(menuname=='_')+3)=[];
        celltypes{cellcountnum}=[menuname ' ' num2str(cellnum)];
        cellcountnum=cellcountnum+1;
    end
end