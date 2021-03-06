load('../granule_traces.mat');
cellnum=14;
ntrials=size(granule_traces,2);
nfake=size(granule_traces,1);
time=length(granule_traces{1,1});

mean_fake=zeros(nfake,length(granule_traces{cellnum,ntrials}));
raw=cell(nfake,1);
for cellnum=1:nfake
    raw{cellnum}=zeros(ntrials,time);
    for i=1:ntrials
        mean_fake(cellnum,:)=mean_fake(cellnum,:)+granule_traces{cellnum,i}/ntrials;
        raw{cellnum}(i,:)=granule_traces{cellnum,i};
    end
end
deadcells=find(max(mean_fake')-min(mean_fake')==0);
mean_fake(deadcells,:)=[];
nfake=size(mean_fake,1);

mean_fake=downsample(mean_fake',10)';
time=floor(time/10)+1;

%% compare dimensionality to real GCs
mean_fake=mean_fake-mean(mean_fake,2)*ones(1,time);
[ef,vf]=eig(mean_fake*mean_fake');
vf=diag(vf);

[~,mean_real]=generate_bases('raw');
nreal=size(mean_real,1);
% mean_real=mean_real(:,1:2001);
mean_real=mean_real(:,500:2500);
mean_real=mean_real-mean(mean_real,2)*ones(1,time);
[er,vr]=eig(mean_real*mean_real');
vr=diag(vr);

%%
h=figure(1);clf;
ndim=8;
subplot(2,1,1);hold on
plot((flipud(vr(end-ndim+1:end)/sum(vr))),'.-')
plot((flipud(vf(end-ndim+1:end)/sum(vf))),'g.-')
legend('real GCs','fake GCs')
title(['Fraction of variance captured by first ' num2str(ndim) ' PCs'])
axis tight
subplot(2,1,2);
semilogy(1:8,(flipud(vr(end-ndim+1:end)/sum(vr))),'.-');
hold on
semilogy(1:8,(flipud(vf(end-ndim+1:end)/sum(vf))),'g.-')
title('Log scaling')
axis tight
fixfigpaper

set(h,'Position',[-700 1 280 300]);
svdir='/Users/Ann/Dropbox/mossyfigs/';
fname = strcat(svdir,'variance_captured_by_PCs');
% print('-depsc','-tiff','-r300',fname);
% saveSameSize(gcf, 'format', 'epsc','file', fname);
%%
h=figure(2);clf;hold on;
for i=1:5
    subplot(5,1,i);
    plot(er(:,end-i+1)'*mean_real);
    hold on
    signflip=sign(mean((er(:,end-i+1)'*mean_real).*(ef(:,end-i+1)'*mean_fake)));
    plot(ef(:,end-i+1)'*mean_fake*signflip,'g');
    axis tight
    if(i==1)
        legend('Real GCs','Fake GCs')
    end
end
subplot(5,1,1);
title('First 5 PCs of responses')
fixfigpaper
set(h,'Position',[-700 1 700 1000]);
fname = strcat(svdir,'first_5_PCs');
print('-depsc','-tiff','-r300',fname);

%% look at fake and real responses sorted along PC1
comp=1;
proj=er(:,end-comp+1)'*mean_real;
fake_sorted=sort_cells(mean_fake,proj);
real_sorted=sort_cells(mean_real,proj);

nplot=60;
offset=1;

h=figure(3);subplot(1,2,1)
imagesc(real_sorted(nplot*(offset-1)+1:nplot*offset,:))
title('Real cells')
subplot(1,2,2);
imagesc(fake_sorted(nplot*(offset-1)+1:nplot*offset,:))
title('Fake cells')
fixfigpaper
set(h,'Position',[-700 1 550 800]);
fname = [svdir 'colormap_of_cell_responses' num2str(offset)];
print('-depsc','-tiff','-r300',fname);

%% plot some responses
% figure(4);clf;hold on
colors='rgbcmkrgbcmk';colors=[colors colors colors];colors=[colors colors colors];
nplot=31;
offset=5;
svdir='/Users/Ann/Dropbox/mossyfigs/';
sc=15;
for offset=1:4
    h=figure(offset);clf;hold on;
    subplot(1,2,1);hold on
    count=1;
    for i=nplot*(offset-1)+1:nplot*offset
        plot(real_sorted(i,:)+i*sc,colors(count));
        count=count+1;
        title('Real cells')
    %     axis([0 2000 0 300]);
        axis tight
    end
    
    subplot(1,2,2);hold on
    count=1;
    for i=nplot*(offset-1)+1:nplot*offset
        plot(fake_sorted(i,:)+i*sc,colors(count));
        count=count+1;
        title('Fake cells')
    %     axis([0 2000 0 300]);
        axis tight
    end
    fixfig
    set(h,'Position',[-1800+400*(offset-1) 1 550 1000]);
    strnum = int2str(offset);
    fname = strcat(svdir,'comparison_of_model_and_granule_cell_subset_number_',strnum);
    print('-depsc','-tiff','-r300',fname);

end


%% test real PCs as a basis for fake cells
p=er^-1*mean_real;
v2=(p*mean_real'*er).^2;
test=er^-1*(mean_fake*mean_fake');
test2=er^-1*(mean_real*mean_real');
figure(1);clf;
subplot(2,1,1);
hold on;
title('Fraction of variance captured by PCs of real cells')
plot(flipud(sum(abs(test2),2)/sum(sum(abs(test2),2))),'.-')
plot(flipud(sum(abs(test),2)/sum(sum(abs(test),2))),'g.-')
legend('Real cells','Fake cells')
xlim([0 50])
subplot(2,1,2)
title('Log scaling')
semilogy(flipud(sum(abs(test2),2)/sum(sum(abs(test2),2))),'.-')
hold on;
semilogy(flipud(sum(abs(test),2)/sum(sum(abs(test),2))),'g.-')
xlim([0 50])
fixfigpaper
set(h,'Position',[-700 1 550 800]);
fname = [svdir 'varcaptured_realPCs'];
print('-depsc','-tiff','-r300',fname);