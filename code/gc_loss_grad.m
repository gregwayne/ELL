function [cst,grad] = gc_loss_grad(theta)

    cst_theta = @(param) gc_loss(param(1),param(2),param(3));
    
    grad = FiniteDifference(cst_theta,theta);
    cst  = cst_theta(theta);

end