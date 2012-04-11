[~,bfns]=generate_bases('thr',.7);
ncells=size(bfns,1);
% bfns=smoothts(bfns,'g',25000,5000);

% [x,y]=meshgrid(linspace(1,4500,ncells),linspace(1,ncells,ncells));
% [xi,yi]=meshgrid(1:4500,1:ncells);
% gbump=normpdf(1:ncells,ncells/2,2);
% bfns=toeplitz(gbump,gbump);
% bfns2=interp2(x,y,bfns,xi,yi);
% clf;imagesc(bfns2)
% bfns=bfns2;
% 

c=bfns*bfns';
% c=eye(ncells);
cinv=(c+1e-5*eye(ncells))^-1;
fbfns=fft(bfns')';

K=fbfns' * cinv * fbfns;


% figure(1);imagesc(linspace(-4500/2,4500/2,4500),linspace(-4500/2,4500/2,4500),fftshift(abs(K))-diag(diag(fftshift(abs(K)))));colorbar
k2=fftshift(real(K))+fftshift(real(K))';

figure(1);imagesc(linspace(-4500/2,4500/2,4500),linspace(-4500/2,4500/2,4500),k2-diag(diag(k2)));colorbar
xlim([-50 50]*4);
ylim([-50 50]*4);
%%
k3=k2-diag(diag(k2));
k4=rot90(k3);
dsig=zeros(4500,1);
dnoise=zeros(4500,1);
ran=-2250:2250;
for i=ran
    dsig(i+2251)=sum(abs(diag(k3,i)));
    dnoise(i+2251)=sum(abs(diag(k4,i)));
end
figure(5);plot(dsig);hold on;plot(dnoise,'g')
%%
figure(4);clf;
plot(ifftshift(diag(k2)));
hold on;
plot(ifftshift(sum(abs(k3),2)),'g')
xlim([0 100])
%%
figure(2);clf;
for rep=1:10
phase=rand(4500/2,1)*2*pi;
amp=zeros(2250,1);
% amp(5)=1;
amp(10)=1;
amp(mod(4500-(1:2250)+1,4500)+1)=amp;
clear phr;
phr=cos(phase)+1i*sin(phase);
phr(mod(4500-(1:2250)+1,4500)+1)=conj(phr);

f = amp.*phr;
% f=zeros(4500,1);
% sptime=ceil(rand(1)*4500);
% f(sptime)=1;
% f=fft(f);
colors='rgbcmkrgbcmkrgbcmk';
figure(2);
subplot(2,1,1);
hold on;plot(abs(K*f),colors(rep));
xlim([0 50])
subplot(2,1,2);
dat=ifft(K*f,'symmetric');
hold on;plot(dat/var(dat),colors(rep));
% plot(sptime,0,'g.')
end