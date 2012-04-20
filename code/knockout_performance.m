[stim resp tran]=loader_plasticity(2);
[stim2 resp2 tran2]=loader_plasticity(1);
[celltypes, bfns]=generate_bases('thr',.8);
rates=convolve_with_synaptic_kernel(bfns',5e-5,0.003,0.02)';
ncells=size(rates,1);
figure(12);clf;hold on;


%CELLTYPES CELLTYPES YES
% 1 = CD_ipsp
% 2 = CDpause_GCs
% 3 = PCA
% 4 = PCA_CDpause_GCs
% 5 = PCA_late
% 6 = PCA_med
% 7 = PCA_med_late

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


rmeans=mean(rates,2);
C=zeros(length(rmeans));
C=rmeans*rmeans';
C=C+1e-5*eye(length(rmeans));
C=C*(4500*mean(mean(rates)))^2;

w=(C^-1*rates*-stim')';

usegrp=cell2mat(cellfun(@(x) x,subpops([noPCA]),'uniformoutput',false));
% usefr=.5;
% usegrp=randperm(ncells);
% usegrp=usegrp(1:round(usefr*ncells));
ratesko=rates;
ratesko(setdiff(1:ncells,usegrp),:)=0;

rmeans=mean(ratesko,2);
C2=zeros(length(rmeans));
C2=rmeans*rmeans';
C2=C2+1e-5*eye(length(rmeans));
C2=C2*(4500*mean(mean(ratesko)))^2;

wko=(C2^-1*ratesko*-stim')'/4;


for i=1:6
    subplot(6,2,i*2-1);
    plot(w(i,:),'k.-');
    hold on;
    plot(find(w(i,:)>0),w(i,w(i,:)>0),'b.');
    plot(find(w(i,:)<0),w(i,w(i,:)<0),'r.');
    xlim([1 124]);
    subplot(6,2,i*2);
    hold on;
    plot(tran,resp(i,:),'g');
    plot(tran,resp2(i,:),'g')
    
    spt=find(stim(i,:));
    plot(tran,w(i,:)*rates,'r');
    plot(tran,wko(i,:)*ratesko,'b');
    plot(tran(find(stim(i,:))),0,'r.')
%     plot(tran,w(i,:)*-ratesko,'c--');
%     plot(tran,wko(i,:)*-rates,'m--');
    xlim([min(tran) max(tran)])
%     ylim([-4 3])

end
%%
figure(15);clf;
colors='rgbcmkrgbcmkrgbcmkrgbcmk';
count=1;
fmax=30;
nf=300;
nph=200;
frange=linspace(1,fmax,nf);
phrange=linspace(0,2*pi,nph);
sig=zeros(nf,nph);
sig2=zeros(nf,nph);
K=C^-1*rates;
for ph=phrange
    count2=1;
    for freq=frange
        stim=cos(2*pi*freq*(1:4500)/4500+ph);
        resp=rates'*(K*stim');
        
%         fresp=abs(fft(resp/sqrt(var(resp))));
%         sig(count2,count)=fresp(round(freq));
        sig2(count2,count)=abs(exp(1i*freq*(1:4500)/4500)*resp(1:end));
        count2=count2+1;
    end
    count=count+1;
end

imagesc(frange,phrange,sig2')

%% compare weights of knockout to original
figure(14);clf;hold on
for i=1:6
    subplot(3,2,i);
    grp=find(wko(i,:));
    plot(w(i,grp),wko(i,grp),'.')
    hold on;
    plot([min(wko(i,grp)) max(wko(i,grp))]*1.25,[min(wko(i,grp)) max(wko(i,grp))]*1.25,'k--')
    xlim([min(wko(i,grp)) max(wko(i,grp))]*1.25)
    ylim([min(wko(i,grp)) max(wko(i,grp))]*1.25)
%     hold on;
%     grp=find(wko(i,:));
%     plot(w(i,grp),'k.-');
%     plot(find(w(i,grp)>0),w(i,grp(w(i,grp)>0)),'b.');
%     plot(find(w(i,grp)<0),w(i,grp(w(i,grp)<0)),'r.');
%     
%     plot(wko(i,grp),'k.-');
%     plot(find(wko(i,grp)>0),wko(i,grp(wko(i,grp)>0)),'c.');
%     plot(find(wko(i,grp)<0),wko(i,grp(wko(i,grp)<0)),'m.');
end

%% sparsify weights to see how sensitive things are
figure(13);clf;hold on;
for i=1:6
    subplot(6,2,i*2-1);
    wsp(i,:)=w(i,:).*(abs(w(i,:))>.2);
    plot(wsp(i,:),'k.-')
    hold on;
    plot(find(wsp(i,:)>0),wsp(i,wsp(i,:)>0),'b.')
    plot(find(wsp(i,:)<0),wsp(i,wsp(i,:)<0),'r.')
    subplot(6,2,i*2);
    hold on;
    plot(tran,resp(i,:),'g');
    plot(tran,resp2(i,:),'g')
    plot(tran,wsp(i,:)*rates);
    plot(tran,w(i,:)*rates,'--');
    plot(tran(find(stim(i,:))),0,'r.')
    xlim([min(tran) max(tran)])
end