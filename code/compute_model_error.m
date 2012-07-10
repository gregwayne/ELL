function [error] = compute_model_error(GC_model,mean_mf,real_cells,errortype)

if(strcmp(errortype,'MSE'))
    [realtrace,modeltrace,tran] = current_based_convolution(GC_model,mean_mf,real_cells);
    realtrace=realtrace-mean(realtrace);
    modeltrace=modeltrace-mean(modeltrace);
    
    error=mean(sqrt((realtrace-modeltrace).^2));
else
    disp('Error: requested errortype undefined.');
end