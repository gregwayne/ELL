[stim resp tran]=loader_plasticity(2);
[stim2 resp2 tran2]=loader_plasticity(1);
[celltypes, bfns]=generate_bases('thr',.8);
rates=convolve_with_synaptic_kernel(bfns',5e-5,0.003,0.02)';
figure(12);clf;hold on;

%CELLTYPES CELLTYPES YES
% 1 = CD_ipsp
% 2 = CDpause_GCs
% 3 = PCA
% 4 = PCA_CDpause_GCs
% 5 = PCA_late
% 6 = PCA_med
% 7 = PCA_med_late

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

w=(C^-1*rates*stim')';

% usegrp=cell2mat(cellfun(@(x) x,subpops([2]),'uniformoutput',false));
% usefr=.5;
% usegrp=randperm(ncells);
% usegrp=usegrp(1:round(usefr*ncells));
ratesko=rates;
ratesko(setdiff(1:ncells,usegrp),:)=0;

rmeans=mean(ratesko,2);
C=zeros(length(rmeans));
C=rmeans*rmeans';
C=C+1e-5*eye(length(rmeans));
C=C*(4500*mean(mean(ratesko)))^2;

wko=(C^-1*ratesko*stim')';


for i=1:6
    subplot(6,2,i*2-1);
    plot(w(i,:),'k.-');
    hold on;
    plot(find(w(i,:)<0),w(i,w(i,:)<0),'r.');
    plot(find(w(i,:)>0),w(i,w(i,:)>0),'b.');
    xlim([1 124]);
    subplot(6,2,i*2);
    hold on;
    plot(tran,resp(i,:),'g');
    plot(tran,resp2(i,:),'g')
    plot(tran,w(i,:)*-rates,'r');
%     plot(tran,wko(i,:)*-ratesko,'b');
%     plot(tran,w(i,:)*-ratesko,'c--');
%     plot(tran,wko(i,:)*-rates,'m--');
    xlim([min(tran) max(tran)])
%     ylim([-4 3])

end

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
    plot(tran,wsp(i,:)*-rates);
    plot(tran,w(i,:)*-rates,'--');
    xlim([min(tran) max(tran)])
end