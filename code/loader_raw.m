function [cellnames, events, recordings] = loader_raw(pth)
% returns a cell array of EOCD times (events) and ephys data (recordings;
% in mV) for all data in folder pth. So, events{1} and recordings{1} would
% be data from the first cell, etc. Time is 5e-5 seconds for both arrays.

folders=dir(pth);
folders = folders(find(cellfun(@isempty,strfind({folders.name},'.'))));
% ^remove things that aren't directories

clear files;
for i=1:length(folders)
    files(i).celltype=folders(i).name;
    files(i).data=dir([pth '/' folders(i).name '/*.mat']);
end

numfiles = sum(cellfun(@length,{files.data}));

count=1;
events=cell(numfiles,1);
recordings=cell(numfiles,1);
for i=1:length(files)
    for j=1:length(files(i).data)
        temp=load([pth '/' folders(i).name '/' files(i).data(j).name]); %load both/all channels into a cell
        datnames = cellfun(@(X) X(end-2:end), fieldnames(temp), 'UniformOutput', false);

        eventind=find(strcmp(datnames,'Ch1'));
        recind=find(strcmp(datnames,'Ch3'));
        backup = find(strcmp(datnames,'Ch4'));

        temp=struct2cell(temp);
        if(~strcmp(temp{recind}.title,'lowgain'))
            recind=backup;
            if(~strcmp(temp{recind}.title,'lowgain'))
                disp('Could not find the input channel!');
            end
        end
        
        name=files(i).celltype;
        name(name=='_')=' ';
        if(strfind(name,'mat'))
            name(strfind(name,'mat')-1:strfind(name,'mat')+2)=[];
        end
        cellnames{count} = [name ' ' num2str(j)];

        events{count}=temp{eventind}.times/temp{recind}.interval + temp{recind}.start;
        recordings{count}=temp{recind}.values;
        if(strcmp(temp{recind}.units,' volt')) %fix units
            recordings{count}=recordings{count}*100;
        end
        count = count+1;
    end
end