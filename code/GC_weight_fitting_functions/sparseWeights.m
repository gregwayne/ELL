function [GC_model,dropCount]=sparseWeights(Wstore,usecell,MSE,w,wtable,weightcost)

GC_model_initialize;

useweights=find(Wstore(usecell,:));
L1                  = sum(abs(w),2);
L0                  = sum(wtable~=0,2);
maxwts=3;
    
%and now compare
% mse_adjusted = mse - mse(1)*((1+bonus)*(ncells-sum(wtable~=0,2)));

mse_adjusted        = MSE + weightcost * L1 + 1e-5 * L0;
mse_adjusted(L1==0) = Inf; %don't return empty fits dude c'mon that helps no one

[~,ind] = min(mse_adjusted);

GC_model.GC_to_model    = usecell;
GC_model.MF_input       = padarray(useweights(nonzeros(wtable(ind,:))),[0 maxwts-sum(wtable(ind,:)~=0)],'post');
GC_model.Ws             = w(ind,:);

dropCount       = sum(wtable(1,:)~=0)-sum(wtable(ind,:)~=0);