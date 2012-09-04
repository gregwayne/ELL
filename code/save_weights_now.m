Wsparse(cellnum,:)=zeros(1,numMFs);
Wsparse(cellnum,nonzeros(GC_model.MF_input)) = nonzeros(GC_model.Ws);

store_sparsities(cellnum) = weightcost;