% [celltypes, bfns]=generate_bases('raw');
% [celltypes, bfns]=generate_bases('2thr',.5, .8); %first is PCA, 2nd is others
[celltypes, bfns]=generate_bases('thr',.7);
ncells=size(bfns,1);

grps=unique(celltypes);
subpops=cell(length(grps),1);
for i=1:length(grps)
	subpops{i}=find(strcmp(celltypes,grps(i)));
end
allgrps=1:7;
mixgrp=4:7;
puregrp=1:3;
nopause=[3 5 6 7];
noPCA=[1 2];
pauseonly=2;

rates=convolve_with_synaptic_kernel(bfns',5e-5,0.003,0.02)';

usegrp=cell2mat(cellfun(@(x) x,subpops([allgrps]),'uniformoutput',false));
rates=rates(usegrp,:);ncells=size(rates,1);



vars=zeros(ncells,1);
NCs=zeros(ncells);


C=zeros(ncells);

[stim, ~, ~]=loader_plasticity(1);
sptime=find(stim(3,:));
win=50;
wmin=max(1,sptime-win);
wmax=min(4500,sptime+win);
C=rates(:,wmin:wmax)*rates(:,wmin:wmax)';

k=makekernel(rates,C);
k=k/sum(rates(:))/4;
% k=k-ones(4500,1)*mean(k,1);



%plot stuff!
tran=(-500:4000)*5e-5;
figure(1);imagesc(tran,1:124,rates);title('basis functions!')
figure(2);imagesc(tran,tran,k);title('kernel!')

check_LR(k,'extraspikes',[]);