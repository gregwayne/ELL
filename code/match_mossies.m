[celltypes,mean_real]=generate_bases('raw');
mean_real=mean_real(:,1:2500);

[mossycelltypes,mossies]=mean_mossy('../final_mossyfibers/');

figure(3);imagesc(smoothts(mossies,'e',200))


%%
ind=find(strcmp(celltypes,'PCA'));
mossyind=find(strcmp(mossycelltypes,'early_mat'));
ind=ind(1:min(30,length(ind)));
mossyind=mossyind(1:min(30,length(mossyind)));

mossysmooth=smoothts([mossies(:,1:200) mossies(:,1:200) mossies(:,1:200) mossies(:,1:200) mossies(:,1:200) mossies],'e',300);
mossysmooth=mossysmooth(:,1001:end);

covcell=(mossysmooth(mossyind,:))*(mean_real(ind,:))';

figure(1);clf;
plot((mean_real(ind,:)-(mean_real(ind,1)-(1:length(ind))'*mean(sqrt(var(mean_real')))*3)*ones(1,2500))');
axis tight
figure(2);clf;
sc=30;
plot((mossysmooth(mossyind,:)-(mossysmooth(mossyind,1)-(1:length(mossyind))'*mean(sqrt(var(mossysmooth')))*3)*ones(1,2500))')
axis tight
% figure(4);imagesc(covcell)
%%
pad=[mossies(:,1:200) mossies(:,1:200) mossies(:,1:200) mossies(:,1:200) mossies(:,1:200)];
mossysmooth=smoothts([pad pad mossies],'e',500);
mossysmooth=mossysmooth(:,2001:end);
figure(3);clf;hold on
plot((mean_real(ind(1),65:end)-mean_real(ind(1),1))/sqrt(var(mean_real(ind(1),:))),'b')
plot((mossysmooth(mossyind(4),:)-mossysmooth(mossyind(4),1))/sqrt(var(mossysmooth(mossyind(4),:)))*1.2,'g')
%%
[celltypes,mean_real]=generate_bases('raw');
mean_real=mean_real(:,1:2500);

thr=.5;
mdf=mean_real(:,2:end)-mean_real(:,1:end-1);
mdf=fft(mdf')';
cutoff=200;
mdf(:,cutoff+2:end-cutoff)=0;
mdf=ifft(mdf')';
mdf=max(mdf,0);
mdf=mdf./(10+max(max(mean_real(:,1:end-1)))-mean_real(:,1:end-1));
mdf=mdf./(sqrt(var(mdf')')*ones(1,length(mdf)))/2;
% mdf(mdf<thr)=0;

spikestore=zeros(size(mean_real));
for cellnum=1:size(mean_real,1)
    [~,i]=findpeaks(mdf(cellnum,:),'minpeakheight',thr);
    figure(1);clf;hold on
    plot(mean_real(cellnum,:)/sqrt(var(mean_real(cellnum,:))));
    plot(mdf(cellnum,:)+3,'g')
    plot(i,mean_real(cellnum,i)/sqrt(var(mean_real(cellnum,:))),'r.')
    pause(.1);
    if(~isempty(i))
        spikestore(cellnum,i)=1;
    end
    spikestore(cellnum,:)=[mdf(cellnum,:) mdf(cellnum,end)];
end
figure(2);imagesc(smoothts(spikestore,'e',50));colormap gray
%%
tran=linspace(-.025,.1,2500);
[~,ind]=max(mossies*smoothts(spikestore,'e',50)');
cellnum=110;

figure(1);clf;hold on;
plot(tran,mean_real(cellnum,:)/sqrt(var(mean_real(cellnum,:))));
plot(tran,spikestore(cellnum,:)/sqrt(var(spikestore(cellnum,:)))/3+2,'r');

ind=find(strcmp(mossycelltypes,'early_mat'));
mnum=ind(37);
plot(tran,smoothts(mossies(mnum,:)/sqrt(var(mossies(mnum,:))),'e',50)+2,'g');
title(num2str(mnum));
xlim([-.005 .021])


%good matches, [GC MF]
% 22    21 or 32
% 4     86 or 91
% 3     84 (kind of)
% 124   73 + 7(ish)
% 110   57 + 37