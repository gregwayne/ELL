cellname = 'pause 9';
[meanid, fid] = find_raw_data(cellname);

mean_trace = struct2cell(load(meanid));
mean_trace = mean_trace{1}.values;

if(isempty(fid))
    break;
end
[events recording] = load_raw_trace(fid);

ncells=length(events);
dat=cell(ncells,1);

minwin=4000;

inter_event_times = events(2:end)-events(1:end-1)-20;
events((inter_event_times<minwin))=[];
inter_event_times(inter_event_times<minwin)=[];

win=floor(min(minwin,min(inter_event_times)));
offset=500;

dat=nan(length(events)-2,win);
for j=2:length(events)-1 %pull out all EOCD events
    stimstart=floor(events(j))-offset;
    temp=recording(stimstart:min(stimstart+offset+win,length(recording)));

    temp(temp>-30)=-30;

    dat(j-1,1:length(temp))=temp;

end

sc=10;

figure(1);clf;
plot(bsxfun(@plus,dat(1:end/sc,:)',(1:size(dat,1)/sc)*10));
title(cellname);
axis tight
figure(2);clf;
plot((-offset:minwin)*5e-5,mean(dat(1:end,:)),'r')
hold on;
tmax=length(mean_trace)-offset;
plot((-offset+1:tmax)*5e-5,mean_trace);

plot([0 0],[min([mean_trace' mean(dat)]) max([mean_trace' mean(dat)])],'k--')
legend('avg from raw traces','avg used for fit')

name=meanid;
ind=strfind(name,'/');
ind=ind(end);
name=name(ind+1:end);
name=regexprep(name,'_',' ');
title(name,'fontsize',15)

xlim([-offset*5e-5 max(minwin,length(mean_trace)-offset)*5e-5]);
fixfig