[events,recordings]=loader_raw('../rawGCdata');
ncells=length(events);
dat=cell(ncells,1);
fdat=cell(ncells,1); %fourier-transformed data (for filtering)

for i=1:ncells %loop over cells
    sptimes{i}=zeros(length(events{i}),20); %store spiketimes (slow)
    win=floor(min(20000,min(events{i}(2:end)-events{i}(1:end-1))-10));
    
    if(i==1||i==8) %events are too close together!
        win=2000;
    end
    dat{i}=nan(length(events{i})-1,win);
    for j=1:length(events{i})-1 %pull out all EOCD events
        stimstart=floor(events{i}(j));
        temp=recordings{i}(stimstart:min(stimstart+win,length(recordings{i})));
        dat{i}(j,1:length(temp))=temp;
        fdat{i}(j,:)=fft(dat{i}(j,:)-mean(dat{i}(j,:)));
        if(~isempty(find(dat{i}(j,:)>0)))
            [pks,locs]=findpeaks(dat{i}(j,:),'minpeakheight',0);
        else
            locs=[];
        end
        sptimes{i}(j,1:length(locs))=locs;
    end
end

%% look to see if it worked
% figure(66);clf;
% ind=5;
% hold on;
% plot(dat{ind}');
% plot(mean(dat{ind}),'k','linewidth',2)
% hold on;
% if(~isempty(nonzeros(sptimes{ind}(:))));
%     plot(nonzeros(sptimes{ind}(:)),0,'g.')
% end