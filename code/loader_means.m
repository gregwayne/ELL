function meanV = loader_means(pth)
%% find all the files
folders=dir(pth);
folders = folders(find(cellfun(@isempty,strfind({folders.name},'.'))));
% ^remove things that aren't directories

clear files;
for i=1:length(folders)-1
    files(i).celltype=folders(i).name;
    files(i).data=dir([pth '/' folders(i).name '/*.mat']);
end

% extract mean responses from files
count=1;
meanV=cell(1,2);
for i=1:length(files)
    for j=1:length(files(i).data)
        temp=struct2cell(load([pth '/' folders(i).name '/' files(i).data(j).name]));
        temp=temp{1};
        name=files(i).celltype;
        name(name=='_')=' ';
        if(strfind(name,'GCs'))
            name(strfind(name,'GCs'):strfind(name,'GCs')+2)=[];
        end
        if(strfind(name,'CD'))
            name(strfind(name,'CD'):strfind(name,'CD')+1)=[];
        end
        name=[name ' ' num2str(j)];
        meanV{count,1}=name;        
        meanV{count,2}=temp.values;
        if(strcmp(temp.units,' volt')) %fix units
            meanV{count,2}=meanV{count,2}*100;
        end
        
        count=count+1;
    end
end