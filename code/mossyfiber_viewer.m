pth='../final_mossyfibers/';
folders=dir(pth);
folders = folders(find(cellfun(@isempty,strfind({folders.name},'.'))));
folders = folders(find(cellfun(@isempty,strfind({folders.name},'modeling'))));

%change these parameters to look at different cells/celltypes.
celltype=2;
cellnum=3;


files=dir([pth folders(celltype).name]);
files(1:2)=[];
files = files(find(cellfun(@isempty,strfind({files.name},'.DS_Store')))); %what why agh

% for cellnum=1:length(files)
data=struct2cell(load(files(cellnum).name));
if(strcmp(data{1}.title,'cmdfilt')||strcmp(data{1}.title,'trigevt')||strcmp(data{1}.title,'cd')) %figure out which channel is which
    eventind=1; spind=2;
elseif(data{1}.length<data{2}.length)
    eventind=1; spind=2;
else
    eventind=2; spind=1;
end

% [data{1}.title ' ' num2str(data{1}.length) ' ' data{2}.title ' ' num2str(data{2}.length) ' ' num2str(eventind)]
figure(1);clf;h=subplot(2,1,1);
hold on;
eventtimes=data{eventind}.times;
spiketimes=data{spind}.times;
tmax=max(eventtimes(2:end)-eventtimes(1:end-1));
tmin=min(eventtimes(2:end)-eventtimes(1:end-1));
meanrsp=zeros(ceil((tmax+.1)/data{eventind}.resolution/100),1);

for i=2:length(eventtimes)-1
    if((eventtimes(i)-eventtimes(i-1)<.2)&&(eventtimes(i+1)-eventtimes(i)>.2)) %only use well-isolated events
        mrkstart=find(spiketimes>eventtimes(i)-.1,1);
        mrkend=find(spiketimes>eventtimes(i)+.3,1); %find the first spike past the end of the interval
        if(mrkstart<mrkend)
            plot((spiketimes(mrkstart:mrkend-1)-eventtimes(i)),i,'b.');
        end
        plot(eventtimes(i+1)-eventtimes(i),i,'g.'); %marks the end of the trial
        spind=ceil((spiketimes(mrkstart:mrkend-1)-eventtimes(i)+.1)/data{eventind}.resolution/100);
        meanrsp(spind)=meanrsp(spind)+1;
    end
end
plot([0 0],[1 length(eventtimes)],'k--')

meanrsp=meanrsp/length(eventtimes)/data{eventind}.resolution/100;
xlim([-.1 .3])
ylim([1 length(eventtimes)-1])
cellname=folders(celltype).name;
cellname(cellname=='_')=' ';
title([cellname ' cell ' num2str(cellnum) ' '])
ylabel('Trial number ')
xlabel('time (s)');
subplot(2,1,2);
hold on
plot(linspace(-.1,tmax,length(meanrsp)),meanrsp,'g');
xlim([-.1 .3])
title('Average response ')
xlabel('time (s) ')
ylabel('firing rate (hz)');
plot(linspace(-.1,tmax,length(meanrsp)),smoothts(meanrsp','g',50,10),'linewidth',2);
legend('raw','smoothed, \sigma^2=2.5ms','Location','northwest')