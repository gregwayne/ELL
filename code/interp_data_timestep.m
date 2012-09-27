function [real_cells, mean_mf] = interp_data_timestep(real_cells,mean_mf,interpfactor)

numGCs=size(real_cells,1);
numMFs=size(mean_mf,1);

%interpolate averaged GC traces
real_temp       = real_cells;
real_cells      = zeros(numGCs,size(real_cells,2)*interpfactor);
for i=1:numGCs
    real_cells(i,:) = interp(real_temp(i,:),interpfactor);
end


%interpolate raw MF trials
% rspstore_temp = rspstore;
% rspstore=cell(numMFs,1);
% for i=1:numMFs
%     rspstore{i} = sparse(zeros(size(mf_temp{i},1),size(mf_temp{i},2)*interpfactor));
% end
% for i=1:numMFs
%     for j=1:size(rspstore{i},1)
%         rspstore{i}(j,:) = sparse(interp(rspstore_temp{i}(j,:),interpfactor));
%     end
% end


%and interpolate the averaged MF traces
mf_temp = mean_mf;
mean_mf = zeros(numMFs,size(mf_temp,2)*interpfactor);
for i=1:numMFs
    mean_mf(i,:) = interp(mf_temp(i,:),interpfactor);
end
