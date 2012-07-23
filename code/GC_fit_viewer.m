%%run this once
load ../GC_fitting_output/GCfitdata_newGCs_fasttau_6inputs.mat; %or what have you

%%
%play with these:
cellstoview=4;
batchnum=1;

tmin=-.025;
tmax=.2;
dt=5e-5;
tran=tmin+dt:dt:tmax;

sparsityfactor=0.4;
dropcount=0;
modeldat_sparse=zeros(1,4500);

posfactor=4;
clc;

for batchnum=1:ceil(numGCs/cellstoview)
fig=figure(1);clf;
for ind=1:cellstoview
    cellnum=(batchnum-1)*cellstoview + ind;
    if(cellnum<=numGCs)
        subplot(cellstoview,2,ind*2-1);
        p=get(gca,'Position');
        p(4)=p(4)*.8;
        set(gca,'Position',p);
        hold on

        [GC_model,dropcount]=sparseWeights(Wstore,cellnum,mean_mf,real_cells,sparsityfactor);
        [~,modeldat_sparse,~] = simulate_current_based_convolution(GC_model,mean_mf,real_cells);
        modeldat_sparse=modeldat_sparse-mean(modeldat_sparse(1:200)); %sparse fits


        modeldat_orig=modeltraces(cellnum,:); %original fits
        modeldat_orig=modeldat_orig-mean(modeldat_orig(1:200));

        realdat=real_cells(cellnum,:)-mean(real_cells(cellnum,1:200));
        
        plot(tran,realdat);
%         plot(tran,modeldat_orig,'r');
        plot(tran,modeldat_sparse,'g');
        xlim([tmin tmax]);
        axis tight
        
        if(ind==3)
            ylabel('Vm relative to pre-EOCD baseline (mV)');
        end
        if(ind==cellstoview)
            xlabel('Time (s)');
            h=legend('Real GC','Model GC (sparse)','Location','NorthWest');
            p=get(h,'Position');
            p(1)=p(1)-.03;
            p(2)=p(2)-.13;
            set(h,'Position',p);
            c=.9;
            set(h, 'Color', [c c c])
            set(h, 'EdgeColor', [c c c])
        end

        tstr=[gctypes{cellnum} '; MSE = ' num2str(MSE(cellnum),'%0.3f')];
        h=title(tstr);
        p=get(h,'Position');
        p(1)=p(1)+.15;
        set(h,'Position',p);

        MFs=find(Wstore(cellnum,:));
        Ws=nonzeros(Wstore(cellnum,:));
        MFs=nonzeros(GC_model.MF_input);
        Ws=nonzeros(GC_model.Ws);
        
        textstr={};
        for i=1:length(MFs)
            pad=blanks(9-length(mftypes{MFs(i)}));
            textstr{i}=[mftypes{MFs(i)} ':' pad 'W = ' num2str(Ws(i),'%0.2f')];
        end
%         textstr{length(textstr)+1}=[num2str(dropcount) ' dropped'];
        ypos=get(gca,'YLim');
        if(abs(ypos(1))>abs(ypos(2)))
            sign=-1;
            pos=min([realdat(2500:4500) modeldat_sparse(2500:4500)]);
        else
            sign=1;
            pos=max([realdat(2500:4500) modeldat_sparse(2500:4500)]);
        end
        if((pos+sign*(ypos(2)-ypos(1))/posfactor)<ypos(1))
            xpos=-.02;
            pos=pos-sign*(ypos(2)-ypos(1))/posfactor;
            ylim([ypos(1)-1.5*(ypos(2)-ypos(1))/posfactor ypos(2)]);
        else
            xpos=0.1;
        end
        h=text(xpos,pos+sign*(ypos(2)-ypos(1))/posfactor,textstr);
        set(h,'FontName','Courier');
        
        subplot(cellstoview,2,ind*2);
        traces=convolve_mossies(GC_model,mean_mf(MFs,:));
        plot(tran,bsxfun(@times,traces',Ws'));
        axis tight;
        p=get(gca,'Position');
        p(4)=p(4)*.8;
        set(gca,'Position',p);
        set(gca, 'YTickLabel', num2str(transpose(get(gca, 'YTick'))))
        box off
        if(ind==cellstoview)
            xlabel('Time (s)');
        end

%         disp(['Inputs fit to ' gctypes{cellnum}]);
%         disp('------');
%         cellfun(@disp,textstr);
%         disp('------');
%         disp(['MSE of fit: ' num2str(MSE(cellnum),'%0.3f')]);
%         disp(' ');
%         disp(' ');
    end
end
set(fig,'PaperPositionMode','manual');
set(fig,'PaperPosition',[0 0 8.5 11]);

fid='../GC_fitting_output/newfits.ps';
if(batchnum==1);
    print(fig, fid, '-dpsc2');
else
    print(fig, fid, '-dpsc2', '-append');
end
end

%%
Wsparse=zeros(size(Wstore));
for cellnum=1:numGCs
    [GC_model,dropcount]=sparseWeights(Wstore,cellnum,mean_mf,real_cells,sparsityfactor);
    Wsparse(cellnum,nonzeros(GC_model.MF_input))=nonzeros(GC_model.Ws);
end
figure;hist(sum(Wsparse~=0,2))
