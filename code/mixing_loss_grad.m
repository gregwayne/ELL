function [cst,grad] = mixing_loss_grad(theta)

    cst_theta = @(param) mixing_loss(param(1),param(2),param(3));
    
    grad = FiniteDifference(cst_theta,theta);
    cst  = cst_theta(theta);

end