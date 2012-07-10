% function [mftypes, rspstore, gctypes, real_cells]=gui_initialize()

[mftypes,rspstore]=export_mossy_rasters('../final_mossyfibers/');
[gctypes,real_cells]=generate_bases('raw');
mean_mf=mean_from_rspstore(rspstore);

%% now let's sort our cells a little bit to make matching easier
% unique_mossies=unique(mftypes_uns);
% 
% for i=1:length(unique_mossies)
%     cellset = find(cell2mat(cellfun(@strcmp,{mftypes_uns},unique_mossies(i),'uniformoutput',false)));
%     marker=zeros(size(cellset));
%     
%     for j=1:length(cellset)
%         meanresponse=mean(rspstore_uns{cellset(j)});
%     
%         if(strcmp(unique_mossies{i},'early'))       %sort early cells by their response onset
%             thr=(min(meanresponse)+max(meanresponse))/2;
%             marker(j)=find(meanresponse>thr,1);
%         elseif(~strcmp(unique_mossies{i},'pause'))  %sort medium and late by the center of their activity
%             thr=(min(meanresponse)+max(meanresponse))/2;
%             marker(j)=mean(find(meanresponse>thr));
%         else                                        %and sort pause cells by their post-pause response onset
%             meanresponse=smoothts(meanresponse,'e',500); %smoothing by a bunch (exponential kernel) lets us find end of the longest inactive interval
%             [~,marker(j)]=min(meanresponse(1000:end));
%         end
%     end
%     
%     [~,ind]=sort(marker);
%     for j=1:length(cellset)
%         rspstore{cellset(j)}=rspstore_uns{cellset(ind(j))};
%         mftypes{cellset(j)}=[mftypes_uns{cellset(ind(j))} ' ' num2str(j)];
%     end
% end
