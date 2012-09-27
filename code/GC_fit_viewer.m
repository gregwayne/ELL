%%run this once
load ../GC_fitting_output/sept10_unrestricted.mat; %or what have you

numGCs=size(Wstore,1);
numMFs=size(Wstore,2);

MSE={};
w={};
wtable={};
cells_to_fit = 1:numGCs;
for usecell = cells_to_fit
    [MSE{usecell}, w{usecell}, wtable{usecell}]    = get_weight_error_tables(Wstore,usecell,mean_mf,real_cells);
end
GC_model_initialize;

%%
%play with these:
cellstoview=4;
tmin    = GC_model.min_t;
tmax    = GC_model.max_t;
dt      = GC_model.dt;
Wfit    = balanced(GC_model);

dropcount = 0;
modeldat_sparse = zeros(1,45000);


savepdf=0;
weightcost=.0005;


for batchnum=1%:ceil(numGCs/cellstoview)
fig=figure(1);clf;
for ind=1:cellstoview
    cellnum=(batchnum-1)*cellstoview + ind;
    if(cellnum<=numGCs)
        subplot(cellstoview,2,ind*2-1);
        p=get(gca,'Position');
        p(4)=p(4)*.8;
        set(gca,'Position',p);
        hold on
        
%         [GC_model,dropcount] = sparseWeights(Wstore,cellnum,MSE{cellnum},w{cellnum},wtable{cellnum},weightcost);
        GC_model = load_weights_from_matrix(GC_model,Wsparse,cellnum);

        MFs=nonzeros(GC_model.MF_input);
        Ws=nonzeros(GC_model.Ws);

        [~,modeldat_sparse,~]   = simulate_current_based_convolution(GC_model,mean_mf,real_cells);
        modeldat_sparse         = modeldat_sparse-mean(modeldat_sparse(1:200)); %sparse fits
        realdat                 = real_cells(cellnum,:)-mean(real_cells(cellnum,1:200));
        
        plot(tran,realdat);
        plot(tran,modeldat_sparse,'g');
        xlim([tmin tmax]);
        axis tight
        
        tstr = [gctypes(cellnum); 'normalized MSE = ' num2str(compute_model_error(GC_model,mean_mf,real_cells,'normMSE'),'%0.3f')];
        h=title(tstr,'interpreter','none');

       
        textstr={};
        possh=4; %usually 4
        for i=1:length(MFs)
            pad=blanks(9-length(mftypes{MFs(i)}));
            textstr{i}=[mftypes{MFs(i)} ':  ' pad num2str(Ws(i)/Wfit,'%0.2f') 'mV'];
        end
        if(savepdf)
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
        end
        
%show the inputs to the model GC
        subplot(cellstoview,2,ind*2);
        traces=convolve_mossies(GC_model,mean_mf(MFs,:));
%         traces=mean_mf(MFs,:);
        plot(tran,zscore(bsxfun(@times,traces',Ws')));
        axis tight;
%         ylim([0 ylimit(2)-ylimit(1)]);
        p=get(gca,'Position');
        p(4)=p(4)*.8;
        set(gca,'Position',p);
        box off
        if(ind==cellstoview)
            xlabel('Time (s)');
        end
        h=title(textstr,'interpreter','none');
        set(h,'fontname','Courier')
    end
end

if(savepdf)
    set(fig,'PaperPositionMode','manual');
    set(fig,'PaperPosition',[0 0 8.5 11]);
    
    fid='../GC_fitting_output/sept10_handsparsed_nocatrestr.ps';
    if(batchnum==1);
        print(fig, fid, '-dpsc2');
    else
        print(fig, fid, '-dpsc2', '-append');
    end
end
end

if 1
    save_weights_now;
end