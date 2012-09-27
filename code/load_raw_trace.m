function [events, recording] = load_raw_trace(pth)


temp=load(pth); %load both/all channels into a cell
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

events=temp{eventind}.times/temp{recind}.interval + temp{recind}.start;
recording=temp{recind}.values;
if(strcmp(temp{recind}.units,' volt')) %fix units
    recording=recording*100;
end