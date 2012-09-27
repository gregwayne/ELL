function [error] = compute_model_error(GC_model,mean_mf,real_cells,errortype)

[realtrace,modeltrace,tran] = simulate_current_based_convolution(GC_model,mean_mf,real_cells);

modeltrace = modeltrace(~isnan(realtrace));
realtrace = realtrace(~isnan(realtrace));

realtrace   = realtrace-mean(realtrace);
modeltrace  = modeltrace-mean(modeltrace);

modeltrace  = modeltrace(~isnan(realtrace));
realtrace   = realtrace(~isnan(realtrace));

if(strcmp(errortype,'MSE'))
    error   = mean((realtrace-modeltrace).^2);
elseif(strcmp(errortype,'normMSE'))
    error   = mean((realtrace-modeltrace).^2)/var(realtrace);
else
    disp('Error: requested errortype undefined.');
end