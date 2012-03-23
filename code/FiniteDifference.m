function grd = FiniteDifference(LossFn,theta)

    eps     = 1e-3;
    grd     = zeros(length(theta),1);
    Id      = eye(length(theta));
    
    for i=1:length(theta)
       
        grd(i) = (LossFn(theta + eps*Id(:,i)) ...
                    - LossFn(theta - eps*Id(:,i)))/(2*eps);
                                
    end
    
        
end