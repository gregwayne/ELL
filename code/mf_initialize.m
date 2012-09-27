function [mftypes,rspstore,mean_mf,MF_indices] = mf_initialize()

[mftypes,rspstore]=loader_mossies('../mossyfibers_final/','activetrials');
fixind = find(1-cellfun(@isempty,strfind(mftypes,'20120419_003')));
rspstore{fixind} = replace_20120419_003('alltrials');

data_dt = 0.05;
mean_mf=zeros(length(rspstore),size(rspstore{1},2)*10);
for i=1:length(rspstore)
    rspstore{i}=rspstore{i}/data_dt;
    mean_mf(i,:) = interp(mean(rspstore{i},1),10);
end


MF_indices = get_struct_of_celltypes(mftypes);
% MF_indices.medlate = MF_indices.med + MF_indices.late;
% MF_indices = rmfield(MF_indices,{'med','late'});

celltypes = fieldnames(MF_indices);
for i=1:length(celltypes) %convert binary vector to indices of cell types
    MF_indices.(celltypes{i})=find(MF_indices.(celltypes{i}));
end