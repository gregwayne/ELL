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

%new pie chart (6/18), with late/UBC as a new category:
nE      = 99;
nM      = 8;
nL      = 2;
nP      = 9;
nEM     = 16;
nEL     = 11;
nEP     = 15;
nML     = 0;
nLP     = 0;
nPM     = 0;
nEMP    = 3;
nELM    = 0;
nELP    = 0;
nLMP    = 0;
nN      = 12;



all_terms_meas      = [nE,nM,nL,nP,nEM,nEL,nEP,nML,nPM,nLP,nELM,nEMP,nELP,nLMP,nN];

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
        for l=0:eps:(1-e-m)
            for p=0:eps:(1-e-m-l)            
                cst = gc_loss_4types(e,m,l,p, all_terms_meas);
                if cst < min_cst

                    min_cst = cst;
                    theta   = [e;m;l;p];

                end

            end
        end
    end
end

[calc] = gc_all_terms_4types(theta(1),theta(2),theta(3),theta(4));
nT          = sum(all_terms_meas);
chi_squared = gc_loss_4types(theta(1),theta(2),theta(3),theta(4), all_terms_meas);
dof         = length(all_terms_meas) - length(theta) - 1;


% Probability that the null hypothesis would generate a less extreme value
% of the chi squared statistic
prob_less   = chi2cdf(chi_squared,dof);
% Probability that the null hypothesis would generate a more extreme value
p_value     = 1-prob_less;

disp(sprintf(strcat('Observation of bins: \n',num2str(all_terms_meas))));
disp(sprintf(strcat('Estimation of bins: \n',num2str(nT*calc))));
disp(sprintf(strcat('Theta vector (e,m,l,p): \n',num2str(theta'))));
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
%%
