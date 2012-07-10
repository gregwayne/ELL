%%run this once
load GCfitdata.mat;

%%
%play with these:
cellstoview=4;
batchnum=1;

tmin=min(tran);
tmax=max(tran);
dt=tran(2)-tran(1);

posfactor=4;
clc;
figure(1);clf;
for i=1:cellstoview
    cellnum=(batchnum-1)*cellstoview + i;
    if(cellnum<=124)
        subplot(cellstoview,1,i);
        hold on
        
        realdat=real_cells(cellnum,:)-mean(real_cells(cellnum,1:200));
        modeldat=modeltraces(cellnum,:)-mean(modeltraces(cellnum,1:200));
        plot(tran,realdat);
        plot(tran,modeldat,'g');
        xlim([tmin tmax]);
        axis tight

        tstr=[gctypes{cellnum} '; MSE = ' num2str(MSE(cellnum),'%0.3f')];
        title(tstr);

        MFs=find(Wstore(cellnum,:));
        Ws=nonzeros(Wstore(cellnum,:));
        textstr={};
        for i=1:length(MFs)
            pad=blanks(9-length(mftypes{MFs(i)}));
            textstr{i}=[mftypes{MFs(i)} ':' pad 'W = ' num2str(Ws(i),'%0.2f')];
        end
        for i=length(MFs)+1:3
            textstr{i}='empty';
        end
        ypos=get(gca,'YLim');
        if(abs(ypos(1))>abs(ypos(2)))
            sign=-1;
            pos=min([realdat(2500:4500) modeldat(2500:4500)]);
        else
            sign=1;
            pos=max([realdat(2500:4500) modeldat(2500:4500)]);
        end
        h=text(.1,pos+sign*(ypos(2)-ypos(1))/posfactor,textstr);
        set(h,'FontName','Courier');
    end
    
    disp(['Inputs fit to ' gctypes{cellnum}]);
    disp('------');
    cellfun(@disp,textstr);
    disp('------');
    disp(['MSE of fit: ' num2str(MSE(cellnum),'%0.3f')]);
    disp(' ');
    disp(' ');
end
fixfig;

%% use this cell to look at the inputs to a given GC
figure(2);clf;

GCnum=6;
MFs=find(Wstore(GCnum,:));
traces=convolve_mossies(GC_model,mean_mf(MFs,:));

%turn this on to view zero-centered traces
% for i=1:length(MFs)
%     traces(i,:)=traces(i,:)-mean(traces(i,1:200));
% end

plot(tran,traces')
xlim([tmin tmax]);
str={};
str{1}=['Input to ' gctypes{GCnum} ' convolved with'];
str{2}= ['\tau_s = ' num2str(GC_model.tau_s) 'ms, \tau_m = ' num2str(GC_model.tau_m) 'ms'];
title(str)
legend(mftypes(MFs))
fixfig