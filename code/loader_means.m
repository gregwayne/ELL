function meanV = loader_means(pth)
%% find all the files

folders=dir(pth);
folders = folders(find(cellfun(@isempty,strfind({folders.name},'.'))));
% ^remove things that aren't directories

clear files;
for i=1:length(folders)
    files(i).celltype=folders(i).name;
    files(i).data=dir([pth '/' folders(i).name '/*.mat']);
end

numfiles = sum(cellfun(@length,{files.data}));

% extract mean responses from files
count=1;
meanV=cell(numfiles,2);
for i=1:length(files)
    for j=1:length(files(i).data)
        expinfo = files(i).data(j).name;
        expdata = struct2cell(load([pth '/' folders(i).name '/' expinfo]));
        expdata = expdata{1};
        
        expinfo = regexp(expinfo,'\((.*)\)','tokens'); %finds the experiment date (which comes between parens in Nate's filenames)
        expinfo = expinfo{:};
        
        name = files(i).celltype;
        name = regexprep(name,'_',' ');
        name = regexprep(name,{'GCs' 'CD'},''); %sometimes there are superfluous strings in directory names
        
        name = [name ' ' num2str(j) ': ' expinfo{:}];
        meanV{count,1} = name;
        meanV{count,2} = expdata.values;
        if(strcmp(expdata.units,' volt')) %fix units
            meanV{count,2} = meanV{count,2}*100;
        end
        
        count=count+1;
    end
end