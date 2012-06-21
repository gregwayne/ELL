function [celltypes meanV]=mean_mossy(pth)
folders=dir(pth);
folders = folders(find(cellfun(@isempty,strfind({folders.name},'.'))));
folders = folders(find(cellfun(@isempty,strfind({folders.name},'modeling'))));

prewin=0.025;

meanV=[];
count=1;
for celltype=1:length(folders)
    files=dir([pth folders(celltype).name]);
    files(1:2)=[];
    files = files(find(cellfun(@isempty,strfind({files.name},'.DS_Store')))); %sometimes this happens
    nsamps=length(files);
    
    meanV=[meanV; zeros(nsamps,2500)];
    for cellnum=1:nsamps
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
        tmax=max(eventtimes(2:end)-eventtimes(1:end-1));
        meanrsp=zeros(ceil((tmax+.4)/5e-5),1); %stored times are in seconds; I convert them to the same scale as the GC data.

        for i=2:length(eventtimes)-1
            if((eventtimes(i)-eventtimes(i-1)>.2)&&(eventtimes(i+1)-eventtimes(i)>.2)) %only use well-isolated events
                mrkstart=find(spiketimes>eventtimes(i)-prewin,1);
                mrkend=find(spiketimes>eventtimes(i+1),1); %find the first spike past the end of the interval
                spind=ceil((spiketimes(mrkstart:mrkend-1)-eventtimes(i)+prewin)/5e-5);
                meanrsp(spind)=meanrsp(spind)+1;
            end
        end
        meanV(count,:)=meanrsp((1:2500))/length(eventtimes);
        celltypes{count}=folders(celltype).name;
        count=count+1;
    end    
end

figure(5);clf;imagesc(smoothts(meanV,'e',50))