function [B,stats] = mylasso(X,Y,varargin)
%LASSO Perform lasso or elastic net regularization for linear regression.
%   [B,STATS] = lasso(X,Y,...) Performs L1-constrained linear least  
%   squares fits (lasso) or L1- and L2-constrained fits (elastic net)
%   relating the predictors in X to the responses in Y. The default is a
%   lasso fit, or constraint on the L1-norm of the coefficients B.
%
%   Positional parameters:
%
%     X                A numeric matrix (dimension, say, NxP)
%     Y                A numeric vector of length N
%   
%   Optional input parameters:  
%
%     'Weights'        Observation weights.  Must be a vector of non-negative
%                      values, of the same length as columns of X.  At least
%                      two values must be positive. (default ones(N,1) or 
%                      equivalently (1/N)*ones(N,1)).
%     'Alpha'          Elastic net mixing value, or the relative balance
%                      between L2 and L1 penalty (default 1, range (0,1]).
%                      Alpha=1 ==> lasso, otherwise elastic net.
%                      Alpha near zero ==> nearly ridge regression.
%     'NumLambda'      The number of lambda values to use, if the parameter
%                      'Lambda' is not supplied (default 100).  Ignored
%                      if 'Lambda' is supplied.  LASSO may return fewer
%                      fits than specified by 'NumLambda' if the residual
%                      error of the fits drops below a threshold percentage 
%                      of the variance of Y.
%     'LambdaRatio'    Ratio between the minimum value and maximum value of
%                      lambda to generate, if the  parameter "Lambda" is not 
%                      supplied.  Legal range is [0,1). Default is 0.0001.
%                      If 'LambdaRatio' is zero, LASSO will generate its
%                      default sequence of lambda values but replace the
%                      smallest value in this sequence with the value zero.
%                      'LambdaRatio' is ignored if 'Lambda' is supplied.
%     'Lambda'         Lambda values. Will be returned in return argument
%                      STATS in ascending order. The default is to have LASSO
%                      generate a sequence of lambda values, based on 'NumLambda'
%                      and 'LambdaRatio'. LASSO will generate a sequence, based
%                      on the values in X and Y, such that the largest LAMBDA                 
%                      value is just sufficient to produce all zero coefficients B.
%                      You may supply a vector of real, non-negative values of 
%                      lambda for LASSO to use, in place of its default sequence.
%                      If you supply a value for 'Lambda', 'NumLambda' and 
%                      'LambdaRatio' are ignored.
%     'DFmax'          Maximum number of non-zero coefficients in the model.
%                      Can be useful with large numbers of predictors.
%                      Results only for lambda values that satisfy this
%                      degree of sparseness will be returned.
%     'Standardize'    Whether to scale X prior to fitting the model
%                      sequence. Note: X and Y are always centered.
%     'RelTol'         Convergence threshold for coordinate descent algorithm.
%                      The coordinate descent iterations will terminate
%                      when the relative change in the size of the
%                      estimated coefficients B drops below this threshold.
%                      Default: 1e-4. Legal range is (0,1).
%     'CV'             If present, indicates the method used to compute MSE.
%                      When 'CV' is a positive integer K, LASSO uses K-fold
%                      cross-validation.  Set 'CV' to a cross-validation 
%                      partition, created using CVPARTITION, to use other
%                      forms of cross-validation. You cannot use a
%                      'Leaveout' partition with LASSO.                
%                      When 'CV' is 'resubstitution', LASSO uses X and Y 
%                      both to fit the model and to estimate the mean 
%                      squared errors, without cross-validation.  
%                      The default is 'resubstitution'.
%     'MCReps'         A positive integer indicating the number of Monte-Carlo
%                      repetitions for cross-validation.  The default value is 1.
%                      If 'CV' is 'resubstitution' or a cvpartition of type
%                      'resubstitution', 'MCReps' must be 1.  If 'CV' is a
%                      cvpartition of type 'holdout', then 'MCReps' must be
%                      greater than one.
%     'PredictorNames' A cell array of names for the predictor variables,
%                      in the order in which they appear in X. 
%                      Default: {}
%     'Options'        A structure that contains options specifying whether to
%                      conduct cross-validation evaluations in parallel, and
%                      options specifying how to use random numbers when computing
%                      cross validation partitions. This argument can be created
%                      by a call to STATSET. CROSSVAL uses the following fields:
%                        'UseParallel'
%                        'UseSubstreams'
%                        'Streams'
%                      For information on these fields see PARALLELSTATS.
%                      NOTE: If supplied, 'Streams' must be of length one.
%   
%   Return values:
%     B                The fitted coefficients for each model. 
%                      B will have dimension PxL, where 
%                      P = size(X,2) is the number of predictors, and
%                      L = length(lambda).
%     STATS            STATS is a struct that contains information about the
%                      sequence of model fits corresponding to the columns
%                      of B. STATS contains the following fields:
%
%       'Intercept'    The intercept term for each model. Dimension 1xL.
%       'Lambda'       The sequence of lambda penalties used, in ascending order. 
%                      Dimension 1xL.
%       'Alpha'        The elastic net mixing value that was used.
%       'DF'           The number of nonzero coefficients in B for each
%                      value of lambda. Dimension 1xL.
%       'MSE'          The mean squared error of the fitted model for each
%                      value of lambda. If cross-validation was performed,
%                      the values for 'MSE' represent Mean Prediction
%                      Squared Error for each value of lambda, as calculated 
%                      by cross-validation. Otherwise, 'MSE' is the mean
%                      sum of squared residuals obtained from the model
%                      with B and STATS.Intercept.
%
%     If cross-validation was performed, STATS also includes the following
%     fields:
%
%       'SE'           The standard error of MSE for each lambda, as
%                      calculated during cross-validation. Dimension 1xL.
%       'LambdaMinMSE' The lambda value with minimum MSE. Scalar.
%       'Lambda1SE'    The largest lambda such that MSE is within 
%                      one standard error of the minimum. Scalar.
%       'IndexMinMSE'  The index of Lambda with value LambdaMinMSE.
%       'Index1SE'     The index of Lambda with value Lambda1SE.
%
%     Examples:
%
%        % (1) Run the lasso on data obtained from the 1985 Auto Imports Database 
%        % of the UCI repository.  
%        % http://archive.ics.uci.edu/ml/machine-learning-databases/autos/imports-85.names
%        load imports-85;
%        Description
%
%        % Extract Price as the response variable and extract non-categorical
%        % variables related to auto construction and performance
%        %
%        X = X(~any(isnan(X(:,1:16)),2),:);
%        Y = X(:,16);
%        Y = log(Y);
%        X = X(:,3:15);
%        predictorNames = {'wheel-base' 'length' 'width' 'height' ...
%            'curb-weight' 'engine-size' 'bore' 'stroke' 'compression-ratio' ...
%            'horsepower' 'peak-rpm' 'city-mpg' 'highway-mpg'};
%
%        % Compute the default sequence of lasso fits.
%        [B,S] = lasso(X,Y,'CV',10,'PredictorNames',predictorNames);
%
%        % Display a trace plot of the lasso fits.
%        axTrace = lassoPlot(B,S);
%        % Display the sequence of cross-validated predictive MSEs.
%        axCV = lassoPlot(B,S,'PlotType','CV');
%        % Look at the kind of fit information returned by lasso.
%        S
%
%        % What variables are in the model corresponding to minimum 
%        % cross-validated MSE, and in the sparsest model within one 
%        % standard error of that minimum.
%        minMSEModel = S.PredictorNames(B(:,S.IndexMinMSE)~=0)
%        sparseModel = S.PredictorNames(B(:,S.Index1SE)~=0)
%
%        % Fit the sparse model and examine residuals.
%        Xplus = [ones(size(X,1),1) X];
%        fitSparse = Xplus * [S.Intercept(S.Index1SE); B(:,S.Index1SE)];
%        corr(fitSparse,Y-fitSparse)
%        figure
%        plot(fitSparse,Y-fitSparse,'o')
%
%        % Consider a slightly richer model. A model with 6 variables may be a 
%        % reasonable alternative.  Find the index for a corresponding fit.
%        df6index = min(find(S.DF==6));
%        fitDF6 = Xplus * [S.Intercept(df6index); B(:,df6index)];
%        corr(fitDF6,Y-fitDF6)
%        plot(fitDF6,Y-fitDF6,'o')         
%         
%        % (2) Run lasso on some random data with 250 predictors
%        %
%        n = 1000; p = 250;
%        X = randn(n,p);
%        beta = randn(p,1); beta0 = randn;
%        Y = beta0 + X*beta + randn(n,1);
%        lambda = 0:.01:.5;
%        [B,S] = lasso(X,Y,'Lambda',lambda);
%        lassoPlot(B,S);
%
%        % compare against OLS
%        %
%        figure
%        bls = [ones(size(X,1),1) X] \ Y;
%        plot(bls,[S.Intercept; B],'.');
%
%        % Run the same lasso fit but restricting the number of
%        % non-zero coefficients in the fitted model.
%        %
%        [B2,S2] = lasso(X,Y,'Lambda',lambda,'DFmax',12);
%
%   See also lassoPlot, ridge, parallelstats.

%   References: 
%   [1] Tibshirani, R. (1996) Regression shrinkage and selection
%       via the lasso. Journal of the Royal Statistical Society,
%       Series B, Vol 58, No. 1, pp. 267-288.
%   [2] Zou, H. and T. Hastie. (2005) Regularization and variable
%       selection via the elastic net. Journal of the Royal Statistical
%       Society, Series B, Vol. 67, No. 2, pp. 301-320.
%   [3] Friedman, J., R. Tibshirani, and T. Hastie. (2010) Regularization
%       paths for generalized linear models via coordinate descent.
%       Journal of Statistical Software, Vol 33, No. 1,
%       http://www.jstatsoft.org/v33/i01.
%   [4] Hastie, T., R. Tibshirani, and J. Friedman. (2008) The Elements
%       of Statistical Learning, 2nd edition, Springer, New York.

%   Copyright 2011 The MathWorks, Inc.

% -------------------------------------- 
% Sanity check the positional parameters
% --------------------------------------

% X a real 2D matrix
if ~ismatrix(X) || length(size(X)) ~= 2 || ~isreal(X)
    error(message('stats:lasso:XnotaReal2DMatrix'));
end

if size(X,1) < 2
    error(message('stats:lasso:TooFewObservations'));
end

% Y a vector, same length as the columns of X
if ~isvector(Y) || ~isreal(Y) || size(X,1) ~= length(Y)
    error(message('stats:lasso:YnotaConformingVector'));
end

% If Y is a row vector, convert to a column vector
if size(Y,1) == 1
    Y = Y';
end

% This screen (okrows) selects all the predictions and response we can use, 
% but there may be further reductions for zero observation weights.
% Defer further checking of X,Y until weights are pre-screened.
okrows = all(isfinite(X),2) & all(isfinite(Y),2);

% --------------------------------------------------------------------
% Parse and process the optional parameters
% --------------------------------------------------------------------

LRdefault = 0.0001;

pnames = { 'weights' 'alpha' 'numlambda' 'lambdaratio' 'lambda' ...
    'dfmax' 'standardize' 'reltol' 'cv' 'mcreps' ...
    'predictornames' 'options' };
dflts  = { []        1       100       LRdefault     []      ...
     []      true          1e-4    'resubstitution'  1 ...
     {}               []};
[weights, alpha, nLambda, lambdaRatio, lambda, ...
    dfmax, standardize, reltol, cvp, mcreps, predictorNames, ParOptions] ...
     = internal.stats.parseArgs(pnames, dflts, varargin{:});

% === 'Alpha' parameter ===

% Require 0 < alpha <= 1.
% "0" would correspond to ridge, "1" is lasso.
if ~isscalar(alpha) || ~isreal(alpha) || ~isfinite(alpha) || ...
        alpha <= 0 || alpha > 1
    error(message('stats:lasso:InvalidAlpha'))
end

% === 'Weights' parameter ===

observationWeights = false;
if ~isempty(weights)
    % This screen works on weights prior to stripping NaNs and Infs.
    if ~isvector(weights) || ~isreal(weights) || size(X,1) ~= length(weights) || ...
            ~all(isfinite(weights)) || any(weights<0) || sum(weights>0) < 2
        error(message('stats:lasso:InvalidObservationWeights'));
    end
    observationWeights = true;
    
    okrows = okrows & weights(:)>0;
    weights = weights(okrows);
    
    % Normalize weights up front.
    weights = weights / sum(weights);
    
    % Below, the convention is that weights is a row vector.
    weights = weights(:)';
end

% Remove observations with NaNs and Infs in the predictor or response
% or with zero observation weight.
X = X(okrows,:);
Y = Y(okrows);

% We need at least two observations after stripping NaNs and Infs.
if size(X,1) < 2
    error(message('stats:lasso:TooFewObservationsAfterNaNs'));
end

% If X has any constant columns, we want to exclude them from the
% coordinate descent calculations.  The corresponding coefficients
% will be returned as zero.
constantPredictors = (range(X)==0);
ever_active = ~constantPredictors;

[N,P] = size(X);

% === 'Standardize' option ===

if ~observationWeights
    muY = mean(Y);
else
    muY = weights*Y;
end
Y0 = bsxfun(@minus,Y,muY);

if standardize
    if ~observationWeights
        % Center and scale
        [X0,~,~] = zscore(X,1);
    else
        % Weighted center and scale
        muX = weights*X;
        X0 = bsxfun(@minus,X,muX);
        sigmaX = sqrt( weights*(X0.^2) );
        % Avoid divide by zero with constant predictors
        sigmaX(constantPredictors) = 1;
        X0 = bsxfun(@rdivide, X0, sigmaX);
    end
else
    if ~observationWeights
        % Center
        muX = mean(X,1);
        X0 = bsxfun(@minus,X,muX);
    else
        % Weighted center
        muX = weights*X;
        X0 = bsxfun(@minus,X,muX);
    end
end

% If using observation weights, make a weighted copy of the 
% predictor matrix, to save time in the weighted partial regressions.

if observationWeights
    wX0 = bsxfun(@times, X0, weights');
end

% Calculate max lambda that permits non-zero coefficients
%
if ~observationWeights
    dotp = abs(X0' * Y0);
    lambdaMax = max(dotp) / (N*alpha);
else
    dotp = abs(sum(bsxfun(@times, wX0, Y0)));
    lambdaMax = max(dotp) / alpha;
end

% === 'Lambda' sequence ===

% Used with nullMSE (calculated below) to terminate
% (ever-less penalized) fits when overfitting is detected.
userSuppliedLambda = true;

if isempty(lambda)
    
    % Used with nullMSE (calculated below) to terminate 
    % (ever-less penalized) fits when overfitting is detected.
    userSuppliedLambda = false;
    
    % Sanity-check of 'NumLambda', should be positive integer.
    if ~isreal(nLambda) || ~isfinite(nLambda) || nLambda < 1
        error(message('stats:lasso:InvalidNumLambda'));
    else
        nLambda = floor(nLambda);
    end
    
    % Sanity-checking of LambdaRatio, should be in [0,1).
    if ~isreal(lambdaRatio) || lambdaRatio <0 || lambdaRatio >= 1
        error(message('stats:lasso:InvalidLambdaRatio'));
    end
    
    if nLambda==1
        lambda = lambdaMax;
    else
        % Fill in a number "nLambda" of smaller values, on a log scale.
        if lambdaRatio==0
                lambdaRatio = LRdefault;
                addZeroLambda = true;
        else
            addZeroLambda = false;
        end
        lambdaMin = lambdaMax * lambdaRatio;
        loghi = log(lambdaMax);
        loglo = log(lambdaMin);
        logrange = loghi - loglo;
        interval = -logrange/(nLambda-1);
        lambda = exp(loghi:interval:loglo)';
        if addZeroLambda
            lambda(end) = 0;
        else
            lambda(end) = lambdaMin;
        end
    end
    
else

    % Sanity check on user-supplied lambda sequence.  Should be non-neg real.
    if ~isreal(lambda) || any(lambda < 0)
        error(message('stats:lasso:NegativeLambda'));
    end

    lambda = sort(lambda(:),1,'descend');
    
end

% === 'RelTol' parameter ===
%
if ~isscalar(reltol) || ~isreal(reltol) || ~isfinite(reltol) || reltol <= 0 || reltol >= 1
    error(message('stats:lasso:InvalidRelTol'));
end

% === 'DFmax' parameter ===
%
% DFmax is #non-zero coefficients 
% DFmax should map to an integer in [1,P] but we truncate if .gt. P
%
if isempty(dfmax)
    dfmax = P;
else
    if ~isscalar(dfmax)
        error(message('stats:lasso:DFmaxBadType'));
    end
    try
        dfmax = uint32(dfmax);
    catch ME
        mm = message('stats:lasso:DFmaxBadType');
        throwAsCaller(MException(mm.Identifier,'%s',getString(mm)));
    end
    if dfmax < 1
        error(message('stats:lasso:DFmaxNotAnIndex'));
    else
        dfmax = min(dfmax,P);
    end
end

% === 'Mcreps' parameter ===
%
if ~isscalar(mcreps) || ~isreal(mcreps) || ~isfinite(mcreps) || mcreps < 1
    error(message('stats:lasso:MCRepsBadType'));
end
mcreps = fix(mcreps);

% === 'CV' parameter ===
%
if isnumeric(cvp) && isscalar(cvp) && (cvp==round(cvp)) && (0<cvp)
    % cvp is a kfold value. Create a cvpartition to pass to crossval. 
    if (cvp > size(X,1))
        error(message('stats:lasso:InvalidCVforX'));
    end
    cvp = cvpartition(size(X,1),'Kfold',cvp);
elseif isa(cvp,'cvpartition')
    if strcmpi(cvp.Type,'resubstitution')
        cvp = 'resubstitution';
    elseif strcmpi(cvp.Type,'leaveout')
        error(message('stats:lasso:InvalidCVtype'));
    elseif strcmpi(cvp.Type,'holdout') && mcreps<=1
        error(message('stats:lasso:InvalidMCReps'));
    end
elseif strncmpi(cvp,'resubstitution',length(cvp))
    % This may have been set as the default, or may have been
    % provided at the command line.  In case it's the latter, we
    % expand abbreviations.
    cvp = 'resubstitution';
else
    error(message('stats:lasso:InvalidCVtype'));
end
if strcmp(cvp,'resubstitution') && mcreps ~= 1
    error(message('stats:lasso:InvalidMCReps'));
end

if isa(cvp,'cvpartition')
    if (cvp.N ~= size(X,1)) || (min(cvp.TrainSize) < 2)
        % We need partitions that match the total number of observations
        % (after stripping NaNs and Infs and zero observation weights), and
        % we need training sets with at least 2 usable observations.
        error(message('stats:lasso:TooFewObservationsForCrossval'));
    end
end

% === 'PredictorNames' parameter ===
%
% If PredictorNames is not supplied or is supplied as empty, we just 
% leave it that way. Otherwise, confirm that it is a cell array of strings.
%
if ~isempty(predictorNames) 
    if ~iscellstr(predictorNames) || length(predictorNames(:)) ~= size(X,2)
        error(message('stats:lasso:InvalidPredictorNames'));
    else
        predictorNames = predictorNames(:)';
    end
end

% === 'Options' parameter ===
% The 'Options' parameter is passed to crossval for handling.
% crossval will do sanity checking.

% --------------------
% Lasso model fits
% --------------------

% Use this to terminate the (ever-less penalized) fits when it becomes
% clear that we are overfitting.
if ~observationWeights
    nullMSE = mean((Y - muY).^2);
else
    % This works because weights are normalized and muY is already weighted
    nullMSE = weights * ((Y - muY).^2);
end

% The struct 'stats' will comprise the second return argument.
% Put place holders for ever-present fields to secure the order
% we want in the struct.
stats = struct();
stats.Intercept      = [];
stats.Lambda         = [];
stats.Alpha          = alpha;
stats.DF             = [];
stats.MSE            = [];
stats.PredictorNames = predictorNames;

[B,Intercept,lambda,mse] = ...
    lassoFit(X,Y, ...
    weights,lambda,alpha,dfmax,standardize,reltol,lambdaMax,ever_active,userSuppliedLambda,nullMSE);

% Store the number of non-zero coefficients for each lambda.
df = sum(B~=0,1);

% ---------------------------------------------------------
% If requested, use cross-validation to calculate 
% Prediction Mean Squared Error for each lambda.
% ---------------------------------------------------------

if ~isequal(cvp,'resubstitution')   
    % Replace dfmax with P, the number of predictors supplied at the
    % command line. dfmax might cause one fold to return empty values, 
    % because no lambda satisfies the dfmax criteria, while other folds 
    % return a numeric value. The lambda sequence has already been 
    % pruned by dfmax, if appropriate, in the call to lassoFit above.
    cvfun = @(Xtrain,Ytrain,Xtest,Ytest) lassoFitAndPredict( ...
        Xtrain,Ytrain,Xtest,Ytest, ...
        lambda,alpha,P,standardize,reltol,lambdaMax,ever_active,true,nullMSE);
    if isempty(weights)
        weights = nan(size(X,1),1);
    end
    cvMSE = crossval(cvfun,[weights(:) X],Y, ...
        'Partition',cvp,'Mcreps',mcreps,'Options',ParOptions);
    mse = mean(cvMSE);
    se  = std(cvMSE) / sqrt(size(cvMSE,1));
    minMSE = min(mse);
    minIx = find(mse==minMSE,1);
    lambdaMin = lambda(minIx);
    minplus1 = mse(minIx) + se(minIx);
    seIx = find((mse(1:minIx) <= minplus1),1,'first');
    if isempty(seIx)
        lambdaSE = [];
    else
        lambdaSE = lambda(seIx);
    end
    
    % Deposit cross-validation results in struct for return value.
    stats.SE           = se;
    stats.LambdaMinMSE = lambdaMin;
    stats.Lambda1SE    = lambdaSE;
    stats.IndexMinMSE  = minIx;
    stats.Index1SE     = seIx;
end

% ------------------------------------------
% Order results by ascending lambda
% ------------------------------------------

nLambda = length(lambda);
reverseIndices = nLambda:-1:1;
lambda = lambda(reverseIndices);
lambda = reshape(lambda,1,nLambda);
B = B(:,reverseIndices);
Intercept = Intercept(reverseIndices);
df = df(reverseIndices);
mse = mse(reverseIndices);
if ~isequal(cvp,'resubstitution')
    stats.SE          = stats.SE(reverseIndices);
    stats.IndexMinMSE = nLambda - stats.IndexMinMSE + 1;
    stats.Index1SE    = nLambda - stats.Index1SE + 1;
end

stats.Intercept = Intercept;
stats.Lambda = lambda;
stats.DF = df;
stats.MSE = mse;
   
end % lasso

% ------------------------------------------
% SUBFUNCTIONS 
% ------------------------------------------

function mse = lassoFitAndPredict(Xtrain,Ytrain,Xtest,Ytest, ...
    lambda,alpha,dfmax,standardize,reltol,lambdaMax,ever_active,userSuppliedLambda,nullMSE)
trainWeights = Xtrain(:,1);
if any(isnan(trainWeights))
    trainWeights = [];
end
Xtrain = Xtrain(:,2:end);
[B,Intercept] = lassoFit(Xtrain,Ytrain, ...
    trainWeights,lambda,alpha,dfmax,standardize,reltol,lambdaMax,ever_active,userSuppliedLambda,nullMSE);
Bplus = [Intercept; B];
testWeights = Xtest(:,1);
if any(isnan(testWeights))
    testWeights = ones(size(Xtest,1),1);
end
Xtest = Xtest(:,2:end);
yfit = [ones(size(Xtest,1),1) Xtest] * Bplus;
mse = testWeights'*(bsxfun(@minus,Ytest,yfit).^2) / sum(testWeights);
end

function [B,Intercept,lambda,mspe] = ...
    lassoFit(X,Y,weights,lambda,alpha,dfmax,standardize,reltol,lambdaMax,ever_active,userSuppliedLambda,nullMSE)
%
% ---------------------------------------------------------------------------------------------------------------------------------
% Perform model fit for each lambda and the given alpha
% ---------------------------------------------------------------------------------------------------------------------------------

[N,P] = size(X);
nLambda = length(lambda);

% If X has any constant columns, we want to exclude them from the
% coordinate descent calculations.  The corresponding coefficients
% will be returned as zero.
constantPredictors = (range(X)==0);
ever_active = ever_active & ~constantPredictors;

% === standardization and weights ===
%
observationWeights = ~isempty(weights);
if ~isempty(weights)
    observationWeights = true;
    weights = weights(:)';
    % Normalize weights up front.
    weights = weights / sum(weights);
end

if ~observationWeights
    muY = mean(Y);
else
    muY = weights*Y;
end
Y0 = bsxfun(@minus,Y,muY);

if standardize
    if ~observationWeights
        % Center and scale
        [X0,muX,sigmaX] = zscore(X,1);
        % Avoid divide by zero with constant predictors
        sigmaX(constantPredictors) = 1;
    else
        % Weighted center and scale
        muX = weights*X;
        X0 = bsxfun(@minus,X,muX);
        sigmaX = sqrt( weights*(X0.^2) );
        % Avoid divide by zero with constant predictors
        sigmaX(constantPredictors) = 1;
        X0 = bsxfun(@rdivide, X0, sigmaX);
    end
else
    if ~observationWeights
        % Center
        muX = mean(X,1);
        X0 = bsxfun(@minus,X,muX);
        sigmaX = 1;
    else
        % Weighted center
        muX = weights*X;
        X0 = bsxfun(@minus,X,muX);
        sigmaX = 1;
    end
end

% If using observation weights, make a weighted copy of the 
% predictor matrix, to save time in the weighted partial regressions.
if observationWeights
    wX0 = bsxfun(@times, X0, weights');
    totalweight = 1;
else
    wX0 = X0;
    totalweight = N;
end

% b will be the current coefficient estimate, iteratively updated.
% Because we retain b from one value of lambda to the next,
% we get a de facto warm start.
b = zeros(P,1);

% Preallocate the returned matrix of coefficients, B, and the intercepts.
B = zeros(P,nLambda);

active = false(1,P);

for i = 1:nLambda
    
    lam = lambda(i);
    if lam >= lambdaMax
        continue;
    end
    threshold = lam * alpha;
    
    % Denominator in coordinate descent update
    if standardize
        if observationWeights
            shrinkFactor = weights*(X0.^2) + lam*(1 - alpha);
        else
            shrinkFactor = ones(1,P) * (1 + lam*(1 - alpha));
        end
    else
        if observationWeights
            shrinkFactor = weights*(X0.^2) + lam*(1 - alpha);
        else
            shrinkFactor = (1/N) * ones(1,N)*(X0.^2) + lam*(1 - alpha);
        end
    end
    
    % Iterative coordinate descent until converged
    while true
        
        bold = b;

        [b,active] = cdescentCycle(X0,wX0,Y0, ...
            b,active,totalweight,shrinkFactor,threshold);
        
        if max( abs(b-bold)./max(1.0,min(abs(b),abs(bold))) ) < reltol
            % Cycling over the active set converged.
            % Do one full pass through the predictors.
            % If there is no predictor added to the active set, we're done.
            % Otherwise, resume the coordinate descent iterations.
            bold = b;
            potentially_active = thresholdScreen(X0,wX0,Y0,b,ever_active,threshold);
            if any(potentially_active)
                new_active = active | potentially_active;
                [b,new_active] = cdescentCycle(X0,wX0,Y0, ...
                    b,new_active,totalweight,shrinkFactor,threshold);
            else
                new_active = active;
            end

            if isequal(new_active, active)
                break
            elseif max( abs(b-bold)./max(1.0,min(abs(b),abs(bold))) ) < reltol
                % We didn't change the coefficients enough by the full pass 
                % through the predictors to warrant continuing the iterations.
                % The active coefficients after this pass and the
                % active coefficients prior to this pass differ.
                % This implies that a coefficient changed between zero
                % and a small numeric value.  There is no "fit-wise" reason
                % to prefer the before and after coefficients, so we
                % choose the more parsimonious alternative.
                if sum(new_active) > sum(active)
                    b = bold;
                else
                    active = new_active;
                end
                break
            else
                active = new_active;
            end
        end
        
    end
    
    B(:,i) = b;
    
    % Halt if maximum model size ('DFmax') has been met or exceeded.
    if sum(active) > dfmax
        % truncate B and lambda output arguments
        lambda = lambda(1:(i-1));
        B = B(:,1:(i-1));
        break
    end
    
    % Halt if we have exceeded a threshold on the percent of
    % residual variance left unexplained.
    if ~userSuppliedLambda
        % Calculate mse of the current fit
        bsig = b ./ sigmaX';
        fit = [ones(size(X,1),1) X] * [(muY-muX*bsig); bsig];
        residuals = bsxfun(@minus, Y, fit);
        if ~observationWeights
            mspe = mean(residuals.^2);
        else
            % This line relies on the weights having been normalized.
            mspe = weights * (residuals.^2);
        end
        if mspe < 1.0e-3 * nullMSE
            lambda = lambda(1:i);
            B = B(:,1:i);
            break
        end
    end
    
end % of lambda sequence

% ------------------------------------------
% Unwind the centering and scaling (if any)
% ------------------------------------------

B = bsxfun(@rdivide, B, sigmaX');
B(~ever_active,:) = 0;
Intercept = muY-muX*B;

% ------------------------------------------
% Calculate Mean Prediction Squared Error
% ------------------------------------------

BwithI = [Intercept; B];
fits = [ones(size(X,1),1) X]*BwithI;
residuals = bsxfun(@minus, Y, fits);
if ~observationWeights
    mspe = mean(residuals.^2);
else
    % This line relies on the weights having been normalized.
    mspe = weights * (residuals.^2);
end

end %-lassoFit

%----------------------------------
% SUBFUNCTIONS
%----------------------------------

% ----------- cdescentCycle() ----------

function [b,active] = cdescentCycle(X0, wX0, Y0, ...
    b, active, totalweight, shrinkFactor, threshold)
%
r = Y0 - X0(:,active)*b(active,:);

for j=find(active);
    bjold = b(j);
    
    % Regress j-th partial residuals on j-th predictor
    rj = r + b(j)*X0(:,j);
    bj = (wX0(:,j)'*rj) / totalweight;
    
    % Soft thresholding
%     b(j) = sign(bj) .* max((abs(bj) - threshold), 0) ./ shrinkFactor(j);                    %soft thresholding performed here!
    b(j) = sign(bj) .* max(((bj) - threshold), 0) ./ shrinkFactor(j); %don't count negative correlations
    if b(j) == 0
        active(j) = false;
    end
    r = r - X0(:,j)*(b(j)-bjold);
end

end %-cdescentCycle

% ===================================================
%                 thresholdScreen() 
% ===================================================

function potentially_active = thresholdScreen(X0, wX0, Y0, ...
    b, active, threshold)
r = Y0 - X0(:,active)*b(active);
% We don't need the (b.*wX2)' term that one might expect, because it
% is zero for the inactive predictors.
potentially_active = abs(r' *wX0) > threshold;                                              %find new cells correlated w/residual
end %-thresholdScreen


