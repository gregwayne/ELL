%%run this once
load ../GC_fitting_output/aug15.mat; %or what have you

numGCs=size(Wstore,1);
numMFs=size(Wstore,2);

%%
cellstoview=4;
batchnum=1;

tmin    = GC_model.min_t;
tmax    = GC_model.max_t;
dt      = GC_model.dt;
tran=tmin+dt:dt:tmax;

sparsityfactor=0.4;
dropcount=0;
modeldat_sparse=zeros(1,4500);

posfactor=4;
clc;

GC_model.V_thresh=-50;
spikes_per_trial=zeros(numGCs,1);

for batchnum=1:ceil(numGCs/cellstoview)
fig=figure(1);clf;
for ind=1:cellstoview
    cellnum=(batchnum-1)*cellstoview + ind;
    if(cellnum<=numGCs)
        subplot(cellstoview,3,ind*3-2);
        p=get(gca,'Position');
        p(4)=p(4)*.8;
        set(gca,'Position',p);
        hold on

%         [GC_model,dropcount]=sparseWeights(Wstore,cellnum,mean_mf,real_cells,sparsityfactor);

%or pull values from Wsparse:
        GC_model.Ws=nonzeros(Wsparse(cellnum,:))';
        GC_model.Ws=padarray(GC_model.Ws,[0 3-length(GC_model.Ws)],'post');
        GC_model.MF_input=find(Wsparse(cellnum,:));
        GC_model.MF_input=padarray(GC_model.MF_input,[0 3-length(GC_model.MF_input)],'post');
        
        [~,modeldat_sparse,~] = simulate_current_based_convolution(GC_model,mean_mf,real_cells);
        modeldat_sparse=modeldat_sparse;%-mean(modeldat_sparse(1:200))+GC_model.E_L; %sparse fits

        realdat=real_cells(cellnum,:)-mean(real_cells(cellnum,1:200))+GC_model.E_L; %not quite right
        
        plot(tran,realdat);hold on;
        plot(tran,modeldat_sparse,'g');
        plot([min(tran) max(tran)],[GC_model.V_thresh GC_model.V_thresh],'k--');

        
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
        

        MFs=find(Wstore(cellnum,:));
        Ws=nonzeros(Wstore(cellnum,:));
        MFs=nonzeros(GC_model.MF_input);
        Ws=nonzeros(GC_model.Ws);
        
        textstr={};
        for i=1:length(MFs)
            pad=blanks(9-length(mftypes{MFs(i)}));
            textstr{i}=[mftypes{MFs(i)} ':' pad 'W = ' num2str(Ws(i),'%0.2f')];
        end

        h=text(0.05,GC_model.V_thresh+5,textstr);
        set(h,'FontName','Courier');
        
        
        xlim([tmin tmax]);
        ylim([-75 -45]);
        
        
        
        
%show the inputs to the model GC
        subplot(cellstoview,3,ind*3-1);
        traces=convolve_mossies(GC_model,mean_mf(MFs,:));
        plot(tran,bsxfun(@times,traces',Ws'));
        axis tight;
        p=get(gca,'Position');
        p(4)=p(4)*.8;
        set(gca,'Position',p);
        set(gca, 'YTickLabel', num2str(transpose(get(gca, 'YTick'))))
        box off
        if(ind==cellstoview)
            h=xlabel('Time (s)');
        end

        %make the raster now so we can talk about it in the title.
        raster{cellnum}=simulate_spike_raster(GC_model,rspstore,real_cells);
        MSE = compute_model_error(GC_model,mean_mf,real_cells,'MSE');
        tstr=[gctypes{cellnum} '; MSE = ' num2str(MSE,'%0.3f') '; avg ' num2str(mean(sum(raster,2))) ' spikes/trial'];
        h=title(tstr);
%         p=get(h,'Position');
%         p(1)=p(1)+.15;
%         set(h,'Position',p);

        spikes_per_trial(cellnum)=mean(sum(raster,2));
        
%show the modelled spike raster
        subplot(cellstoview,3,ind*3);hold on;
        p=get(gca,'Position');
        p(4)=p(4)*.8;
        set(gca,'Position',p);
        box off
        for i=1:size(raster,1);
            if(~isempty(find(raster{cellnum}(i,:))))
                plot(tran(find(raster{cellnum}(i,:))),i,'b.');
            end
        end
        xlim([GC_model.min_t*1e-3, GC_model.max_t*1e-3]);
        ylim([1 100]);
        
        
%and some raw traces
%         subplot(cellstoview,4,ind*4);hold on;
%         p=get(gca,'Position');
%         p(4)=p(4)*.8;
%         set(gca,'Position',p);
%         for i=1:4
%             [~,modeldat,~]=simulate_current_based_expeuler(GC_model,rspstore,real_cells);
%             plot(tran,modeldat,colors(i));
%         end
%         xlim([-.025 .2]);
%         ylim([-65 GC_model.V_thresh+10]);
        
        
    end
end

set(fig,'PaperPositionMode','manual');
set(fig,'PaperPosition',[0 0 8.5 11]);

fid='../GC_fitting_output/spiking_m50.ps';
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
