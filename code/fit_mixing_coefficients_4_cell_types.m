function theta = fit_mixing_coefficients_4_cell_types(Wsparse,mftypes)
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
% nE      = 99;
% nM      = 8;
% nL      = 2;
% nP      = 9;
% nEM     = 16;
% nEL     = 11;
% nEP     = 15;
% nML     = 0;
% nLP     = 0;
% nPM     = 0;
% nEMP    = 3;
% nELM    = 0;
% nELP    = 0;
% nLMP    = 0;
% nN      = 12;


all_terms_meas      = get_mixing_input_counts(Wsparse,mftypes,4);
all_terms_meas(1) = all_terms_meas(1)+31;
all_terms_meas(5) = all_terms_meas(5)+1;

nN = 21; %what should this be?
all_terms_meas = [all_terms_meas nN];

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
eps     = 0.025;
min_p   = 1e-4;

min_cst = inf;

for e=min_p:eps:1
    for m=min_p:eps:(1-e)
        for l=min_p:eps:(1-e-m)
            for p=min_p:eps:(1-e-m-l)            
                cst = mixing_loss_4types(e,m,l,p, all_terms_meas);
                if cst < min_cst                    
                    min_cst = cst;
                    theta.e=e;
                    theta.m=m;
                    theta.l=l;
                    theta.p=p;
                end

            end
        end
    end
end

[calc]      = mixing_all_terms_4types(theta.e,theta.m,theta.l,theta.p);
nT          = sum(all_terms_meas);
chi_squared = mixing_loss_4types(theta.e,theta.m,theta.l,theta.p, all_terms_meas);
dof         = length(all_terms_meas) - length(fieldnames(theta)) - 1;

% Probability that the null hypothesis would generate a less extreme value
% of the chi squared statistic
prob_less   = chi2cdf(chi_squared,dof);
% Probability that the null hypothesis would generate a more extreme value
p_value     = 1-prob_less;

% disp(sprintf(strcat('Observation of bins: \n',num2str(all_terms_meas))));
% disp(sprintf(strcat('Estimation of bins: \n',num2str(nT*calc))));
% disp(sprintf(strcat('Theta vector (e,m,l,p): \n',num2str(theta'))));
% disp(sprintf('Degrees of Freedom: %d',dof));
% disp(sprintf('Chi-Squared Statistic: %f',chi_squared));
% disp(sprintf('p-value: %f',p_value));
% 
% if p_value < 0.05
%     disp(sprintf('Reject Null Hypothesis: Random Mixing Model Matches Poorly'));
% else
%     disp(sprintf('Retain Null Hypothesis: Random Mixing Model Matches Well Enough'));
% end


% bar_matrix = [nT * calc; all_terms_meas];                                                                                                       
% bar(bar_matrix);                                                                                                                      
% bar(bar_matrix');                                                                                                                     
% set(gca,'XTickLabel',{'Early', 'Medium', 'Late', 'Pause', 'E+M', 'E+L', 'E+P', 'M+L', 'P+M', 'L+P', 'E+L+M', 'E+M+P', 'E+L+P', 'L+M+P','None'})
% legend('Fit', 'Measured')
