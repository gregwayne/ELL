close all;

% Grid Search to Optimize Chi-Squared Statistic
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

[calc,meas] = gc_all_terms(theta(1),theta(2),theta(3));
nT          = sum(meas);
chi_squared = gc_loss(theta(1),theta(2),theta(3));
dof         = length(meas) - length(theta) - 1;

% Probability that the null hypothesis would generate a less extreme value
% of the chi squared statistic
prob_less   = chi2cdf(chi_squared,dof);
% Probability that the null hypothesis would generate a more extreme value
p_value     = 1-prob_less;

disp(sprintf(strcat('Observation of bins: \n',num2str(meas),'\n')));
disp(sprintf(strcat('Estimation of bins: \n',num2str(nT*calc),'\n')));
disp(sprintf(strcat('Theta vector (e,m,p): \n',num2str(theta'),'\n')));
disp(sprintf('Degrees of Freedom: %d \n',dof));
disp(sprintf('Chi-Squared Statistic: %f \n',chi_squared));
disp(sprintf('p-value: %f \n',p_value));

if p_value < 0.05
    disp(sprintf('Reject Null Hypothesis: Random Mixing Model Matches Poorly'));
else
    disp(sprintf('Retain Null Hypothesis: Random Mixing Model Matches Well Enough'));
end