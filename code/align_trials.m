[cellnames, events,recordings]=loader_raw('../raw_gcs_mat');
ncells=length(events);
dat=cell(ncells,1);

minwin=4000;

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
        
    end
end

%% look to see if it worked

cellnum=38;
gcexp = regexp(gctypes{cellnum},': ([0-9_]*[0-9_])','tokens');
ind = find(1-cellfun(@isempty,strfind(cellnames,gcexp{:}{:})));
ind=ind(1);

figure(33);clf;
% hold on;
% imagesc(dat{ind}(1:end-1,:));

sc=1;
plot(bsxfun(@plus,dat{ind}(1:sc:end-1,:)',(1:ceil(size(dat{ind},1)/sc-1))*5));
title(cellnames(ind),'interpreter','none');
axis tight
figure(2);clf;
plot(mean(dat{ind}(1:end-1,:)),'k','linewidth',2)
hold on;
plot((1:length(real_cells(cellnum,:)))*.1,real_cells(cellnum,:));