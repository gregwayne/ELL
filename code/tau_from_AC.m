%run align_trials to get raw data
align_trials;
%%

count=1;colors='rgbcmkrgbcmkrgbcmk';
figure(1);clf;hold on;
figure(2);clf;hold on;

usecells=[5 6 9 10];
ncells=length(usecells);
paramfit=zeros(ncells,3);
paramvar=zeros(ncells,3);
paramstore=[];
convwin=50;
tran=0:GC_model.dt:GC_model.max_t;

%exponential function we're fitting:
fh = @(x,p) p(1) + p(2)*exp(-x./p(3));
%error function:
errfh = @(p,x,y) sum((y(:)-fh(x(:),p)).^2);
% an initial guess of the exponential parameters:
p0 = [3 50 1];

for cellnum=usecells
    sample=dat{cellnum};
    ntrials=size(sample,1);
    bump=1;
    cutoff=size(sample,2)-bump;
    convstore=zeros(ntrials,cutoff-5);
    sampleplot=zeros(ntrials,cutoff);
    
    params=zeros(ntrials,3);
    for trial=1:ntrials
        test=sample(trial,bump+1:bump+cutoff);
        test=test-mean(test);

        % drop spikes
        if(~isempty(find(test>5)))
            win=150;
            [~,ind]=findpeaks(test,'minpeakheight',5,'minpeakdistance',win/2);
            for i=1:length(ind)
                test(max(ind(i)-win,1):min(ind(i)+win,length(test)))=test(max(ind(i)-win-1,1));
                test(test<-10)=0;
            end
        end
        sampleplot(trial,:)=test;

        dttest=test(3:end)-test(1:end-2);
        dttest=dttest-mean(dttest);

        dtconv=conv(dttest,fliplr(dttest));
        dtconv(1:ceil(length(dtconv)/2)+2)=[];
        dtconv=dtconv*sign(mean(dtconv(1:5)));
        convstore(trial,:)=dtconv;
        
        %fit an exponential
        rel=(convstore(trial,2:convwin)+convstore(trial,1:convwin-1))/2;
        t=tran(1:convwin-1);
        params(trial,:) = fminsearch(errfh,p0,[],t,rel);
    end
    paramstore=[paramstore; params];
    paramfit(count,:)=mean(params);
    paramvar(count,:)=var(params);
    
    figure(1);
    subplot(ceil(length(usecells)/2),2,count);
    hold on;
    plot(((convstore(:,2:convwin)+convstore(:,1:convwin-1)))'/2)
    expfit=paramfit(count,1)+paramfit(count,2)*exp(-tran(1:convwin)/paramfit(count,3));
    plot(expfit,'k','linewidth',2)
    axis tight
    figure(2);hold on;
    plot(mean(convstore(:,2:convwin)+convstore(:,1:convwin-1))'/2/mean2(convstore(:,1:2)),colors(count))
    axis tight

    count=count+1;
end
paramfit

%%
figure(1);clf;
subplot(2,1,1);cla;hold on;
plot(test-mean(test));  %plot Vm
plot(smoothts(dttest*10,'g',20,4),'g')     %plot time derivative
axis tight
subplot(2,1,2);cla;hold on;
plot(dtconv,'.-')
hold on;
plot(expfit,'g','Linewidth',2);
