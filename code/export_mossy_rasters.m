function [celltypes,rspstore]=export_mossy_rasters(pth)
% pth='../final_mossyfibers/';
folders=dir(pth);
folders = folders(find(cellfun(@isempty,strfind({folders.name},'.'))));
folders = folders(find(cellfun(@isempty,strfind({folders.name},'modeling'))));

winstart=-.025;
winend=.2;
dt=5e-5;

rspstore=cell(94,1); %hacks!

cellcountnum=1;
for celltype=1:4
    files=dir([pth folders(celltype).name]);
    files(1:2)=[];
    files = files(find(cellfun(@isempty,strfind({files.name},'.DS_Store')))); %what why agh
    for cellnum=1:length(files)
        data=struct2cell(load([pth folders(celltype).name '/' files(cellnum).name]));
        if(strcmp(data{1}.title,'cmdfilt')||strcmp(data{1}.title,'trigevt')||strcmp(data{1}.title,'cd')) %figure out which channel is which
            eventind=1; spind=2;
        elseif(data{1}.length<data{2}.length)
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
        celltypes{cellcountnum}=folders(celltype).name;
        cellcountnum=cellcountnum+1;
    end
end