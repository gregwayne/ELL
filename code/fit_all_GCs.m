data_initialize();
GC_model_initialize();

numGCs=size(real_cells,1);
numMFs=size(mean_mf,1);

%%

Wstore=zeros(numGCs,numMFs);
modeltraces=zeros(numGCs,(GC_model.max_t-GC_model.min_t)/GC_model.dt);
MSE=zeros(numGCs,1);

%set some parameters for our two fitting functions
DFmax=6;

%these 3 are just if we're doing a grid search:
wbounds=[0 1000];
gridres=4;
nreps=8;

clear stats;

for modelcell=1:numGCs %loop over all the real GC's!
    disp(['Fitting cell ' num2str(modelcell)]);
    GC_model.GC_to_model = modelcell;
    
    % pick mossy inputs!
    [GC_model,~] = fitterlasso(GC_model,mean_mf,real_cells,DFmax);

    ncells = length(nonzeros(GC_model.MF_input));
    wmin = wbounds(1)*ones(ncells,1);
    wmax = wbounds(2)*ones(ncells,1);

    % fit weights!
%     GC_model = fitterLS(GC_model,mean_mf,real_cells);
% 	GC_model = fittergrid(GC_model,mean_mf,real_cells,wmin,wmax,wbounds,gridres,nreps);
    

    % and save our results
    [~,modeltraces(modelcell,:),~] = simulate_current_based_convolution(GC_model,mean_mf,real_cells);
    MSE(modelcell) = compute_model_error(GC_model,mean_mf,real_cells,'MSE');
    w=zeros(1,numMFs);
    w(nonzeros(GC_model.MF_input)) = GC_model.Ws(find(GC_model.MF_input));
    Wstore(modelcell,:)= w;
end
tran=GC_model.min_t+GC_model.dt:GC_model.dt:GC_model.max_t;
%%
save('GCfitdata_newGCs_fasttau_5inputs','GC_model','mean_mf','real_cells','tran','modeltraces','MSE','Wstore','gctypes','mftypes');