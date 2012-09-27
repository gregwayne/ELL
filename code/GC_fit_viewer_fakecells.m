%%run this once, or load the data from a file
GCparams = load_GC_generation_parameters();
[mftypes,rspstore,mean_mf,MF_indices] = mf_initialize();

% load ~/Dropbox/GC_sim_params;

%%
cellstoview = 4;
cellstomodel = 200;
batchran = 1:ceil(cellstomodel/cellstoview);

GC_model_initialize;
tmin    = GC_model.min_t;
tmax    = GC_model.max_t;
dt      = GC_model.dt;
tran    = tmin+dt:dt:tmax;

posfactor=4;
savepdf = 1;

for batchnum=batchran
fig=figure(1);clf;
    for ind=1:cellstoview
        
        GC_model = make_fake_GC_2(GCparams,MF_indices);
        
        subplot(cellstoview,3,ind*3-2);
        p=get(gca,'Position'); p(4)=p(4)*.8; set(gca,'Position',p);
        hold on
        if(ind==3)
            ylabel('Vm relative to pre-EOCD baseline (mV)');
        end
        if(ind==cellstoview)
            xlabel('Time (s)');
        end
        
        [~,modeldat_sparse,~] = simulate_current_based_convolution(GC_model,mean_mf,[]);
        plot(tran,modeldat_sparse);
        hold on;
        plot([min(tran) max(tran)],[GC_model.V_thresh GC_model.V_thresh],'k--');

        xlim([tmin tmax]);
        ylim([-75 -35]);




    %show the inputs to the model GC
        subplot(cellstoview,3,ind*3-1);
        p=get(gca,'Position'); p(4)=p(4)*.8; set(gca,'Position',p);
        set(gca, 'YTickLabel', num2str(transpose(get(gca, 'YTick'))))
        box off
        if(ind==cellstoview)
            xlabel('Time (s)');
        end
        textstr={};
        for i=1:3
            if(GC_model.MF_input(i))
                textstr{i}=[mftypes{GC_model.MF_input(i)} ':    W = ' num2str(GC_model.Ws(i),'%0.2f')];
            end
        end
        title(textstr,'interpreter','none');

        traces=convolve_mossies(GC_model,mean_mf(nonzeros(GC_model.MF_input),:));
        plot(tran,bsxfun(@times,traces',nonzeros(GC_model.Ws)'));
        axis tight;
        drawnow;






        %make the raster now so we can talk about it in the title.
        raster=simulate_spike_raster(GC_model,rspstore,50);

    %show the modelled spike raster
        subplot(cellstoview,3,ind*3);hold on;
        p=get(gca,'Position');
        p(4)=p(4)*.8;
        set(gca,'Position',p);
        box off
        for i=1:size(raster,1);
            if(~isempty(find(raster(i,:))))
                plot(tran(find(raster(i,:))),i,'b.');
            end
        end
        xlim([GC_model.min_t, GC_model.max_t]);
        ylim([1 50]);

    end

    set(fig,'PaperPositionMode','manual');
    set(fig,'PaperPosition',[0 0 8.5 11]);

    if(savepdf)
        fid='../GC_fitting_output/spiking_new.ps';
        if(batchnum==1);
            print(fig, fid, '-dpsc2');
        else
            print(fig, fid, '-dpsc2', '-append');
        end

    end
end

%%
Wsparse=zeros(size(Wstore));
for cellnum=1:numGCs
    [GC_model,dropcount]=sparseWeights(Wstore,cellnum,mean_mf,real_cells,sparsityfactor);
    Wsparse(cellnum,nonzeros(GC_model.MF_input))=nonzeros(GC_model.Ws);
end
figure;hist(sum(Wsparse~=0,2))
