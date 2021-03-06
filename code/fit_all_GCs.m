%load data
GC_model_initialize;
data_initialize();


%% fit inputs with lasso

Wstore              = zeros(numGCs,numMFs);
% modeltraces         = zeros(numGCs,(GC_model.max_t-GC_model.min_t)/GC_model.dt);
% MSE                 = zeros(numGCs,1);
% tran                = GC_model.min_t+GC_model.dt:GC_model.dt:GC_model.max_t;
DFmax               = 10; %max number of cells lasso is allowed to fit

% if matlabpool('size')==0
%     matlabpool('open');
% end

GC_models={};
for i=1:numGCs
    GC_models{i}    = GC_model;
end
MF_indices          = get_struct_of_celltypes(mftypes);
GC_indices          = get_struct_of_celltypes(gctypes);

cells_to_fit        = 1:numGCs;
tic
for modelcell=cells_to_fit
    disp(['Fitting cell ' num2str(modelcell)]);

    GC_model = GC_models{modelcell};
    GC_model.GC_to_model = modelcell;
    
    inputs_to_use = regexp(gctypes{modelcell},'[a-z]*[a-z]?[ a-z]*[a-z]','match');
    inputs_to_use = regexp(inputs_to_use{:},'\s','split');
    indarray = MF_indices.(inputs_to_use{1});
    for i=2:length(inputs_to_use);
        indarray = indarray + MF_indices.(inputs_to_use{i});
    end
    
%     indarray = ones(numMFs,1);
    
    % pick mossy inputs
    [GC_model,~] = fitterlasso(GC_model,mean_mf,real_cells,indarray,DFmax);

    % fit weights
%     GC_model = fitterLS(GC_model,mean_mf,real_cells);
    

    % and save our results
%     [~,modeltraces(modelcell,:),~] = simulate_current_based_convolution(GC_model,mean_mf,real_cells);
%     MSE(modelcell) = compute_model_error(GC_model,mean_mf,real_cells,'MSE');

    if(size(GC_model.MF_input,1)~=1)
        GC_model.MF_input=GC_model.MF_input';
    end
    if(size(GC_model.Ws,1)~=1)
        GC_model.Ws=GC_model.Ws';
    end
    GC_models{modelcell} = GC_model;
end
toc

for modelcell = cells_to_fit
    ind=nonzeros(GC_models{modelcell}.MF_input.*(GC_models{modelcell}.Ws>0));
    Wstore(modelcell,ind) = nonzeros(GC_models{modelcell}.Ws);
end

%% sparsify the fit weight matrix Wstore to make Wsparse
% old code: doesn't work with our new sparseWeights implementation (for tweaking sparseness costs)

% Wfit=balanced(GC_model);
% weightcost=1/(10*Wfit); %weight on L1 penalty (there's also an L0 penalty but we have that weight (=1e-5) hard-coded in right now)
% Wsparse=zeros(size(Wstore));
% tic
% for cellnum=cells_to_fit
%     cellnum
%     [GC_model,dropcount]=sparseWeights(Wstore,cellnum,mean_mf,real_cells,weightcost);
%     [nonzeros(GC_model.MF_input); nonzeros(GC_model.Ws)]
%     Wsparse(cellnum,nonzeros(GC_model.MF_input))=nonzeros(GC_model.Ws);
% end
% toc
% figure;hist(sum(Wsparse~=0,2))
%%
fid='../GC_fitting_output/sept10_unrestricted';
save(fid,'GC_model','mean_mf','real_cells','Wstore','gctypes','mftypes','MSE','w','wtable','store_sparsities','Wsparse');