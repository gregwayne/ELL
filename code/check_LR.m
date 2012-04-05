function check_LR(k)
%compare kernel performance to observed negative image effect

[stim resp tran]=loader_plasticity();
ndat=size(stim,1);

z=zeros(ndat,4500);
figure(3);clf;
hold on;
for i=1:ndat
    subplot(ceil(ndat/2),2,i);
    plot(tran,stim(i,:));
    hold on;
    plot(tran,resp(i,:),'g');
    z(i,:)=-k*stim(i,:)';
    z(i,:)=z(i,:)/sqrt(var(z(i,:)));
    
    z(i,:)=z(i,:)-mean(z(i,:)); %scaling/centering of z: kosher?
    
    plot(tran,z(i,:),'r');
    xlim([0 max(tran)])
    ylim([-4 3])
end
subplot(ceil(ndat/2),2,8);
imagesc(tran,tran,k)

% subplot(ceil(ndat/2),2,8);
figure(4);clf;
hold on;
flip=[-1 -1 1 1 1 1 1]';
for i=1:ndat
    plot(tran,flip(i)*(resp(i,:)-z(i,:)),'color',[0 i/ndat 1-i/ndat]);
end
plot(tran,mean((flip*ones(1,4500)).*(resp-z)),'k--','linewidth',2)
xlim([0 max(tran)]);
title('kernel - observed')