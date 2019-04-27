%% calcValuations
% Monte-Carlo simulation of payoffs under geometric brownian motion for
% external assets. Uses the Eisenberg/Noe algorithm to compute payoffs.
%
% Eisenberg, L., & Noe, T. H. (2001). Systemic risk in financial systems. 
% Management Science, 47(2), 236-249.
%
% *Inputs*
%
% * vecE: vector (banks x 1) of external assets
% * vecMu: vector (banks x 1) of mean returns of external assets
% * matCovariance: variance-covariance matrix (banks x banks) of banks'
% external assets
% * matL: matrix (banks x banks) of interbank claims
% * numSimulations: number of simulation paths for Monte-Carlo
% * dT (optional): time increment. Default 1
% * strCalculationMethod (optional): 'Standard' for Eisenberg/Noe
% algorithm (default), 'Iterate' for iterative solution
%
% *Outputs*
%
% * matValuations: (banks x numSimulations) matrix of payoff vectors under
% different simulation paths. Payoffs are not discounted yet!
%


function matValuations = calcValuations(vecE,vecMu,matCovariance,matL,numSimulations,dT,strCalculationMethod)

%% Get inputs & Declarations
if nargin < 7 ; strCalculationMethod = 'Standard'; end
if nargin < 6 ; dT = 1; end

numBanks = length(vecE);
vecVariance = diag(matCovariance);

%% Simulate external asset returns
L = chol(matCovariance,'lower');
dW = L*randn(numBanks,numSimulations);
matE = repmat(vecE, 1, numSimulations);
matE = matE.*exp((vecMu-0.5*vecVariance)*dT+sqrt(dT)*dW);

%% Compute payoffs for each simulation path
matValuations = nan(numBanks,numSimulations);

for i = 1:numSimulations
    vecE_Final = matE(:,i);
    vecPayments = calcPayments(vecE_Final,matL,strCalculationMethod);
    matValuations(:,i) = vecPayments;
end

end