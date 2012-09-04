[cellnames, events,recordings]=loader_raw('../raw_gcs_mat');
ncells=length(events);
dat=cell(ncells,1);

minwin=4000;

% fdat=cell(ncells,1); %fourier-transformed data (for filtering)

for i=1:ncells %loop over cells
%     sptimes{i}=zeros(length(events{i}),20); %store EOCD times (slow)

    inter_event_times = events{i}(2:end)-events{i}(1:end-1)-20;
    events{i}((inter_event_times<2000))=[];
    inter_event_times(inter_event_times<minwin)=[];
    
    win=floor(min(minwin,min(inter_event_times)));
    
    dat{i}=nan(length(events{i})-1,win);
    for j=2:length(events{i})-1 %pull out all EOCD events
        stimstart=floor(events{i}(j))-500;
        temp=recordings{i}(stimstart:min(stimstart+500+win,length(recordings{i})));
        dat{i}(j-1,1:length(temp))=temp;
        
%         fdat{i}(j,:)=fft(dat{i}(j,:)-mean(dat{i}(j,:)));
%         if(~isempty(find(dat{i}(j,:)>0)))
%             [pks,locs]=findpeaks(dat{i}(j,:),'minpeakheight',0);
%         else
%             locs=[];
%         end
%         sptimes{i}(j,1:length(locs))=locs;
    end
end

%% look to see if it worked
ind = 1;
figure(1);clf;
% hold on;
imagesc(dat{ind}(1:end-1,:));
%     plot(mean(dat{ind}),'k','linewidth',2)
% 
% xlim([0 4500])
% ylim([-90 0])
title(cellnames(ind));

