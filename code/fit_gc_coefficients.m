addpath minFunc/
options.display = 'on';
options.method  = 'lbfgs';
options.maxIter = 1000;

theta0          = theta;

[theta_opt,cst] = minFunc(@(par) gc_loss_grad(par),theta0);

eps     = 0.01;

min_cst = inf;

for e=0:eps:1
    for m=0:eps:(1-e)
        for p=0:eps:(1-e-m)            

            cst = gc_loss(e,m,p);
            if cst < min_cst
                
                min_cst = cst;
                theta   = [e;m;p];
                
            end
            
        end
    end
end

theta
min_cst