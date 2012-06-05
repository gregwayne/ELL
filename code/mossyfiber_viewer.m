pth='../mossyfibers/';
folders=dir(pth);
folders = folders(find(cellfun(@isempty,strfind({folders.name},'.'))));

%%
%change these parameters to look at different cells/celltypes.
celltype=1;
cellnum=4;

files=dir([pth folders(celltype).name]);
files(1:2)=[];
data=struct2cell(load(files(cellnum).name));
if(strcmp(data{1}.title,'cmdfilt')) %figure out which channel is which
    eventind=1; spind=2;
else
    eventind=2; spind=1;
end



figure(1);clf;hold on;
eventtimes=round(data{eventind}.times/data{eventind}.resolution);
spiketimes=round(data{spind}.times/data{spind}.resolution);
mrk=find(spiketimes>eventtimes(1),1);
nsp=zeros(length(eventtimes)-1,1);
for i=1:length(eventtimes)-1
    mrkend=find(spiketimes>eventtimes(i+1),1); %find the first spike past the end of the interval
    if(mrk<mrkend)
        plot((spiketimes(mrk:mrkend-1)-eventtimes(i))*data{spind}.resolution,i,'b.');
        nsp(i)=mrkend-mrk;
    end
    mrk=mrkend;
end
ylim([1 length(eventtimes)-1])
cellname=folders(celltype).name;
cellname(cellname=='_')=' ';
title([cellname ' cell ' num2str(cellnum)])
ylabel('Trial number')
xlabel('is this in seconds?')
% fixfig %if you're fancy

figure(2);hist(nsp,40)
title(['Average of ' num2str(mean(nsp)) ' spikes per trial']);
% fixfig