function [cst,grd] = FitCost(VS,PS,theta)

    f       = @(theta) L(VS,PS,theta);    
    cst     = f(theta);
    cst
    grd     = FiniteDifference(f,theta);
    
end