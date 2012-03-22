function [meanV] = loader_means(pth)
%% find all the files
folders=dir(pth);
folders = folders(find(cellfun(@isempty,strfind({folders.name},'.'))));
% ^remove things that aren't directories

clear files;
for i=1:length(folders)-1
    files(i).celltype=folders(i).name;
%     cd(folders(i).name);
    files(i).data=dir([pth '/' folders(i).name '/*.mat']);
%     cd ..;
end

% extract mean responses from files
count=1;
meanV=cell(1);
for i=1:length(files)
%     cd([pth '/' files(i).celltype]);
    for j=1:length(files(i).data)
        temp=struct2cell(load([pth '/' files(i).celltype '/' files(i).data(j).name]));
        temp=temp{1};
        meanV{count}=temp.values;
        count=count+1;
    end
end