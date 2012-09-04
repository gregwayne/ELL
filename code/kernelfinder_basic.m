% pth='../DATA2';
% rcell=loader_means(pth);
tmax=2000;

rmat=cell2mat(cellfun(@(x) [x(1:min(length(x),tmax))' x(end)*ones(1,max(tmax-length(x),0))],rcell(:,2),'UniformOutput',false));

%%
rmat=convolve_mossies(GC_model,bfn(find(sum(bfn,2)),:));
rmat=rmat(:,1:tmax);

%%

rmat=rmat-mean(rmat,2)*ones(1,tmax);
ncells=size(rmat,1);
C=rmat*rmat'/tmax^2 + diag(ones(ncells,1))*1e-5;

%let's make up some noise correlations!
% nc_same=0;
% nc_diff=0;
% nc_self=0;
% NC=zeros(ncells);
% for i=1:ncells
%     NC(i,i)=nc_self;
%     for j=i+1:ncells
%         if(strcmp(rcell(i,1),rcell(j,1)))
%             NC(i,j)=nc_same; NC(j,i)=nc_same;
%         else
%             NC(i,j)=nc_diff; NC(j,i)=nc_diff;
%         end 
%     end
% end
% 
% C=C+NC;
Cinv=C^-1;
K=rmat'*Cinv*rmat;
D=zeros(tmax*2-1,1);
for i=1:tmax*2-1
    D(i)=sum(diag(K,i-tmax));
end

% let's take a look!
figure(1);
subplot(2,1,1);
imagesc(K)
subplot(2,1,2);
plot(-tmax+1:tmax-1,D);
%%
fig=figure(1);clf;
for ind=1:6
    subplot(6,2,ind*2-1)
    plot(tran,resp(ind,:));
    hold on
    plot(tran,resp(ind,:)*K/2e7,'g');
    xlim([min(tran) max(tran)]);
    if(ind==1)
        title('Fit to negative image');
    end
    
    subplot(6,2,ind*2);
    plot(tran,-stim(ind,:));
    hold on;
    plot(tran,-stim(ind,:)*K/2e5,'g');
    xlim([min(tran) max(tran)]);
    if(ind==1)
        title('Fit to pulse');
    end
end
set(fig,'PaperPositionMode','manual');
set(fig,'PaperPosition',[0 0 8.5 11]);
fid='../GC_fitting_output/fitplasticity.ps';
print(fig, fid, '-dpsc2');

%%
fig=figure(2);clf;
plot(tran,bsxfun(@plus,convolve_mossies(GC_model,bfn)',(size(bfn,1):-1:1)/1000))
title('Averaged spiking activity of model GCs')
axis tight
set(fig,'PaperPositionMode','manual');
set(fig,'PaperPosition',[0 0 8.5 11]);
fid='../GC_fitting_output/bfns_traces.ps';
print(fig, fid, '-dpsc2');

fig=figure(3);clf;
imagesc(tran,1:size(bfn,1),convolve_mossies(GC_model,bfn))
title('Averaged spiking activity of model GCs')
set(fig,'PaperPositionMode','manual');
set(fig,'PaperPosition',[0 0 8.5 11]);
fid='../GC_fitting_output/bfns_image.ps';
print(fig, fid, '-dpsc2');