%% calcValuations
% Monte-Carlo simulation of payoffs under geometric brownian motion for
% external assets. Uses the Rogers/Veraart algorithm to compute payoffs.
%
% Rogers, L. C., & Veraart, L. A. (2013). Failure and rescue in an
% interbank network. Management Science, 59(4), 882-898.
%
% *Inputs*
%
% * vecE: vector (banks x 1) of external assets
% * vecMu: vector (banks x 1) of mean returns of external assets
% * matCovariance: variance-covariance matrix (banks x banks) of banks'
% external assets
% * matL: matrix (banks x banks) of interbank claims
% * numSimulations: number of simulation paths for Monte-Carlo
% * dblAlpha (optional): recovery value on other assets (default 1)
% * dblBeta (optional): recovery value on interbank assets (default 1)
% * dT (optional): time increment. Default 1
%
% *Outputs*
%
% * matValuations: (banks x numSimulations) matrix of payoff vectors under
% different simulation paths. Payoffs are not discounted yet!
%


function matValuations = calcValuations(vecE,vecMu,matCovariance,matL,numSimulations,dblAlpha,dblBeta,dT)

%% Get inputs & Declarations
if nargin < 8 ; dT = 1; end

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
    vecPayments = calcPayments(vecE_Final,matL,dblAlpha,dblBeta);
    matValuations(:,i) = vecPayments;
end

end