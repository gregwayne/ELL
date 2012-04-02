%fake respose that uses the low-frequency components of the mean signal,
%plus random noise with the same power spectrum as the high-frequency terms
%of non-averaged data. Cutoff determines where the low-freq/high-freq
%boundary lies (in arbitrary not-hertz units because recordings are odd
%lengths). Pretty dumb, because power spectrum of the response changes with
%time.

% align_trials(); %call this once to load dat
cellnum=7;
test=fft(mean(dat{cellnum}));
test2=zeros(1,size(dat{cellnum},2));
for i=1:size(dat{cellnum},1)
    test2=test2+abs(fft(dat{cellnum}(i,:)-mean(dat{cellnum}(i,:))))/size(dat{cellnum},1);
end
cutoff=20;
test(cutoff:end-cutoff+2)=0;


phrnd=rand(1,floor(length(cutoff+1:length(test)-cutoff+1)/2))*2*pi;
if(mod(length(test),2)==0)
    phrnd=[fliplr(phrnd) 0 -(phrnd)];
else
    phrnd=[fliplr(phrnd) -(phrnd)];
end

noise=test2(cutoff+1:end-cutoff+1).*(cos(phrnd)+1i*sin(phrnd));
noise=[zeros(1,cutoff) noise zeros(1,cutoff-1)];

% for comparison: raw data, mean response (green), fake response (blue)
clf;hold on;
plot(dat{cellnum}');
plot(ifft(test+noise),'b','linewidth',2)
plot(mean(dat{cellnum}),'g','linewidth',2)
