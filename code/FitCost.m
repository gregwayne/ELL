function [cst,grd] = FitCost(VS,PS,theta)

    LossFn  = @(par) L(VS,PS,par);    
    cst     = LossFn(theta);
    grd     = FiniteDifference(LossFn,theta);
    
end