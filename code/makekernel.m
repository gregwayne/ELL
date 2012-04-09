function [K] = makekernel(rates,C)
% rates = firing rates of cells (NxT)
% C = covariance matrix (cov(rates) + individual var + NCs)


if(det(C)~=0)
    Cinv=C^-1;
else
    eps=1e-5;
    Cinv = (C+eye(length(C))*eps)^-1;
%     disp('warning, Cinv was close to singular');
end

K = rates'*Cinv*rates;