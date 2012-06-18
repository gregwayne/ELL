close all;

%%
%old pie chart
% nE      = 101;
% nM      = 4;
% nP      = 7;
% nEM     = 20;
% nEP     = 12;
% nPM     = 0;
% nEMP    = 0;
% nN      = 14;

%new pie chart (6/18)
nE      = 110;
nM      = 8;
nP      = 9;
nEM     = 16;
nEP     = 15;
nPM     = 0;
nEMP    = 3;
nN      = 14; %12 from nate's email + two late/UBC alone cells

%if we were to make late/UBC a new category:
% nE      = 99;
% nM      = 8;
% nL      = 2;
% nP      = 9;
% nEM     = 16;
% nEL     = 11; %seems high given the number of late-only cells
% nEP     = 15;
% nML     = 0;
% nLP     = 0;
% nPM     = 0;
% nEMP    = 3;
% nELM    = 0;
% nELP    = 0;
% nLMP    = 0;
% nELMP   = 0;
% nN      = 12;



all_terms_meas      = [nE,nM,nP,nEM,nEP,nPM,nEMP,nN];
%    s = sum(all_terms_meas);
% 
%    all_terms_meas = zeros(1, 8);
%    for i = 1:s
%        idx = randi(8);
%        all_terms_meas(idx) = all_terms_meas(idx) + 1;
%    end
% 
%    disp(all_terms_meas);



% Grid Search to Optimize Chi-Squared Statistic
eps     = 0.01;

min_cst = inf;

for e=0:eps:1
    for m=0:eps:(1-e)
        for p=0:eps:(1-e-m)            

            cst = gc_loss(e,m,p, all_terms_meas);
            if cst < min_cst
                
                min_cst = cst;
                theta   = [e;m;p];
                
            end
            
        end
    end
end

[calc] = gc_all_terms(theta(1),theta(2),theta(3));
nT          = sum(all_terms_meas);
chi_squared = gc_loss(theta(1),theta(2),theta(3), all_terms_meas);
dof         = length(all_terms_meas) - length(theta) - 1;

% Probability that the null hypothesis would generate a less extreme value
% of the chi squared statistic
prob_less   = chi2cdf(chi_squared,dof);
% Probability that the null hypothesis would generate a more extreme value
p_value     = 1-prob_less;

disp(sprintf(strcat('Observation of bins: \n',num2str(all_terms_meas))));
disp(sprintf(strcat('Estimation of bins: \n',num2str(nT*calc))));
disp(sprintf(strcat('Theta vector (e,m,p): \n',num2str(theta'))));
disp(sprintf('Degrees of Freedom: %d',dof));
disp(sprintf('Chi-Squared Statistic: %f',chi_squared));
disp(sprintf('p-value: %f',p_value));

if p_value < 0.05
    disp(sprintf('Reject Null Hypothesis: Random Mixing Model Matches Poorly'));
else
    disp(sprintf('Retain Null Hypothesis: Random Mixing Model Matches Well Enough'));
end


bar_matrix = [nT * calc; all_terms_meas];                                                                                                       
bar(bar_matrix);                                                                                                                      
bar(bar_matrix');                                                                                                                     
set(gca,'XTickLabel',{'Early', 'Medium', 'Pause', 'Early + Medium', 'Early + Pause', 'Pause + Medium', 'E + M + P', 'None'})                              
legend('Fit', 'Measured')
