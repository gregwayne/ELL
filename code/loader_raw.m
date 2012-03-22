function [events, recordings] = loader_raw(pth)
% returns a cell array of EOCD times (events) and ephys data (recordings;
% in mV) for all data in folder pth. So, events{1} and recordings{1} would
% be data from the first cell, etc. Time is 5e-5 seconds for both arrays.

files=dir([pth '/*.mat']); %load all mat files
events=cell(length(files),1);
recordings=cell(length(files),1);
for i=1:length(files)
    temp=struct2cell(load([pth '/' files(i).name])); %load both channels into a cell
    if(strcmp(temp{1}.title,'lowgain')) %figure out which channel is which
        stimind=1; eventind=2;
    else
        stimind=2; eventind=1;
    end
    
    events{i}=temp{eventind}.times/temp{stimind}.interval + temp{stimind}.start;
    recordings{i}=temp{stimind}.values;
    if(strcmp(temp{stimind}.units,' volt')) %fix units
        recordings{i}=recordings{i}*100;
    end
end