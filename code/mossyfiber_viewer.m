pth='../mossyfibers/';
folders=dir(pth);
folders = folders(find(cellfun(@isempty,strfind({folders.name},'.'))));
folders = folders(find(cellfun(@isempty,strfind({folders.name},'modeling'))));
%%
%change these parameters to look at different cells/celltypes.
celltype=1;
cellnum=1;
% for celltype=1:4
files=dir([pth folders(celltype).name]);
files(1:2)=[];
files = files(find(cellfun(@isempty,strfind({files.name},'.DS_Store')))); %what why agh

% for cellnum=1:length(files)
data=struct2cell(load(files(cellnum).name));
if(strcmp(data{1}.title,'cmdfilt')) %figure out which channel is which
    eventind=1; spind=2;
else
    eventind=2; spind=1;
end



figure(1);clf;h=subplot(2,1,1);
hold on;
eventtimes=data{eventind}.times;
spiketimes=data{spind}.times;
mrk=find(spiketimes>eventtimes(1),1);
nsp=zeros(length(eventtimes)-1,1);
tmax=max(eventtimes(2:end)-eventtimes(1:end-1));
tmin=min(eventtimes(2:end)-eventtimes(1:end-1));
meanrsp=zeros(ceil((tmax+.1)/data{eventind}.resolution/100),1);

ran=10.^(-10:.001:0);
n=zeros(size(ran));

for i=1:length(eventtimes)-1
    mrkstart=find(spiketimes>eventtimes(i)-.1,1);
    mrkend=find(spiketimes>eventtimes(i+1),1); %find the first spike past the end of the interval
    if(mrkstart<mrkend)
        plot((spiketimes(mrkstart:mrkend-1)-eventtimes(i)),i,'b.');
        n=n+hist(spiketimes(mrkstart+1:mrkend)-spiketimes(mrkstart:mrkend-1),ran);
        nsp(i)=mrkend-mrk;
    end
    plot(eventtimes(i+1)-eventtimes(i),i,'g.'); %marks the end of the trial
    spind=ceil((spiketimes(mrkstart:mrkend-1)-eventtimes(i)+.1)/data{eventind}.resolution/100);
    meanrsp(spind)=meanrsp(spind)+1;
    mrk=mrkend;
end
plot([0 0],[1 length(eventtimes)],'k--')

meanrsp=meanrsp/length(eventtimes)/data{eventind}.resolution/100;
xlim([-.1 tmin])
ylim([1 length(eventtimes)-1])
cellname=folders(celltype).name;
cellname(cellname=='_')=' ';
title([cellname ' cell ' num2str(cellnum) ' '])
ylabel('Trial number ')
xlabel('time (s)');
subplot(2,1,2);
hold on
plot(linspace(-.1,tmax,length(meanrsp)),meanrsp,'g');
xlim([-.1 tmin])
title('Average response ')
xlabel('time (s) ')
ylabel('firing rate (hz)');
fixfigpaper
plot(linspace(-.1,tmax,length(meanrsp)),smoothts(meanrsp','g',50,10),'linewidth',2);
legend('raw','smoothed, \sigma^2=2.5ms','Location','northwest')


% svdir='/Users/Ann/Dropbox/mossyfigs/';
% fname=[cellname ' cell ' num2str(cellnum)];
% fname(fname==' ')='_';
% print('-depsc','-tiff','-r300',[svdir fname]);
% 
% end
% end

%%
ran=10.^(-10:.01:0);
n=hist(spiketimes(2:end)-spiketimes(1:end-1),ran);
figure(2);clf;bar(log(ran),log(n))
xlim([-7 -1])
%%
set(gcf,'units','pixels','Position',[0 0 300 400]);
set(gcf,'paperpositionmode','auto')

svdir='/Users/Ann/Dropbox/mossyfigs/';
fname=[folders(celltype).name '_cell' num2str(cellnum)];
saveSameSize(gcf, 'format', 'png','file', [svdir fname]);