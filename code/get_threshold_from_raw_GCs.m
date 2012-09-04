function [threshold_mean,threshold_var] = get_threshold_from_raw_GCs(pth)

%load raw GC traces
[~,~,raw_gcs] = loader_raw(pth);
numGCs=length(raw_gcs);
spikes=cell(numGCs,1);

for i=1:numGCs
    [spikes{i},~]=find_spike_indices(raw_gcs{i});
end

%% get threshold from raw traces

abs_thr=zeros(numGCs,1);
baseline=zeros(numGCs,1);
for i=1:numGCs
    if(~isempty(spikes{i})) %only look at GCs that spiked
    
        %first compute the spike-triggered average membrane potential
        win=300;
        STA=zeros(win*2+1,1);
        for j=1:length(spikes{i})
            STA=STA+raw_gcs{i}(spikes{i}(j)-win:spikes{i}(j)+win)/length(spikes{i});
        end
    
        %we tried several methods of finding the threshold (using first or
        %second derivatives), but the waveform of spikes varied so much
        %that no one method worked universally. So we settled on just
        %taking the value of the membrane potential 15 timesteps (0.75ms)
        %prior to the peak of the action potential wave form-- which seemed
        %to give reasonable values, by eye.
        [~,ind]=max(STA);
        back=15;
        abs_thr(i)=STA(ind-back);
        
        %this method of determining the baseline could also stand
        %improvement, as it's sensitive to extreme perturbations in the
        %recording.
        baseline(i)=min(raw_gcs{i});
    end
end

rel_thr=abs_thr-baseline;
rel_thr(rel_thr>50)=0;      %some of our values seem extreme so I'm tossing them out
rel_thr=nonzeros(rel_thr);

threshold_mean = mean(rel_thr);
threshold_var  = var(rel_thr);