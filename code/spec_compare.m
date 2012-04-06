%compare the average power spectrum to the power spectrum of the average,
%to see what frequency components lose power when we average. Try
%(hamfistedly) to correct for effect of spikes on the power spectrum.

cellnum=5;
meanspec=abs(fft(mean(dat{cellnum})-mean(mean(dat{cellnum}))));

trspec=zeros(size(meanspec));
for i=1:size(dat{cellnum},1)
    trspec=trspec+abs(fft(dat{cellnum}(i,:)-mean(dat{cellnum}(i,:))))/size(dat{cellnum},1);
end

trfix=trspec;
spikeslope=(mean(trfix(300:305))-mean(trfix(200:205)))/100;
yint = mean(trfix(300:305))-spikeslope*300;
trsub=(1:floor(length(trfix)/2))*spikeslope + yint;
trsub=([trsub(1) trsub fliplr(trsub)]);
trfix=trfix-trsub(1:length(trfix))*.9;
trfix(300:length(trfix)-300+2)=0;

figure(2);
clf;
subplot(2,1,1)
plot(dat{cellnum}')
subplot(2,1,2);
plot(meanspec);hold on;
plot(trspec,'g');
plot(trfix,'r');
xlim([0 200])
ylim([-10 3000])

