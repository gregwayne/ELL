function [stim resp tran]=loader_plasticity(day)
switch day
    case 1
        pth='../dataformodel/plasticity/20111213';
    otherwise
        pth='../dataformodel/plasticity/20111220';
end
files = dir(pth);
files = files(find(1-cellfun(@isempty,strfind({files.name},'.mat'))));
nfiles=length(files);

tmax=4500;
stim=zeros(nfiles,tmax);
resp=zeros(nfiles,tmax);
stimtimes=zeros(nfiles,1);
for i=1:nfiles
    temp=struct2cell(load([pth '/' files(i).name]));
    temp=temp{1};
    stim(i,:)=0;
    stimtime=round(max((temp.delay-temp.start)/temp.interval,1));
    stim(i,stimtime)=100; %should it be temp.delay-temp.start or temp.delay?
    stimtimes(i)=temp.delay-temp.start;
    resp(i,:)=temp.values;
end
tstep=temp.interval;

ind=[stimtimes (1:nfiles)'];
ind=sortrows(ind,1);ind=ind(:,2);
stim=stim(ind,:);
resp=resp(ind,:);
tran=(1:tmax)*tstep;