% function check_LR(k,varargin)
%compare kernel performance to observed negative image effect

[stim resp tran]=loader_plasticity(2);
[stim2 resp2 tran2]=loader_plasticity(1);
ndat=size(stim,1);
dt=tran(2)-tran(1);

z=zeros(ndat,4500);
figure(3);clf;
hold on;
for i=1:ndat
%     if(~isempty(varargin))
%         if(strcmp(varargin{1},'extraspikes'))
%             for j=1:length(varargin{2})
%                 stim(i,round((varargin{2}(j)-tran(1))/dt)+1)=1;
%             end
%         end
%     end
    subplot(ceil(ndat/2),2,i);
    plot(tran,stim(i,:));
    hold on;
    plot(tran,resp(i,:),'g');
    if(i<size(stim2,1))
        plot(tran,resp2(i,:),'g');
    end
    
%windowed correlations!%%%%%%%%%%%
    [celltypes, bfns]=generate_bases('thr',.7);
    rates=convolve_with_synaptic_kernel(bfns',5e-5,0.001,0.01)';
    sptime=find(stim(i,:));
    win=10;
    wmin=max(1,sptime-win);
    wmax=min(4500,sptime+win);
    C=rates(:,wmin:wmax)*rates(:,wmin:wmax)';
%     C=sum(rates,2)*sum(rates,2)';
%     C=zeros(size(C));
    k=makekernel(rates,C);%/sum(rates(:))/3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    z(i,:)=-k*stim(i,:)';
    z(i,:)=z(i,:)/sqrt(var(z(i,:)));

    z(i,:)=z(i,:)-mean(z(i,:)); %scaling/centering of z: kosher? decidedly not.
    
    plot(tran,smooth(z(i,:),300),'r');
    xlim([min(tran) max(tran)])
    ylim([-4 3])
end
% subplot(ceil(ndat/2),2,8);
% imagesc(tran,tran,k)

% subplot(ceil(ndat/2),2,8);
% figure(4);clf;
% hold on;
% for i=1:ndat
%     plot(tran,resp(i,:)-z(i,:),'color',[0 i/ndat 1-i/ndat]);
% end
% plot(tran,mean(resp-z),'k--','linewidth',2)
% xlim([0 max(tran)]);
% title('kernel - observed')