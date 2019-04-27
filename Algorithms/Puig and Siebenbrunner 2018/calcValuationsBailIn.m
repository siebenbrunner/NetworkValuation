%% calcValuationsBailIn
% Monte-Carlo simulation of payoffs under geometric brownian motion for
% external assets. Uses the Puig/Siebenbrunner algorithm with bail-in
% to compute payoffs.
%
%
% *Inputs*
%
% * vecE: vector (banks x 1) of other assets
% * vecMu: vector (banks x 1) of mean returns of external assets
% * matCovariance: variance-covariance matrix (banks x banks) of banks'
% external assets
% * matL: matrix (banks x banks x seniorities) of interbank claims
% * matTheta: matrix (banks x banks) of interbank holdings
% * numK: integer, determines the number of seniority classes corresponding to bail-in 
% * funConversion: function that returns the desired conversion factor.
%               Inputs: MatBailIn - matrix (banks x numK) of bail-in amount per bank and seniority class
%                       vecEquity: equity vector (banks x 1)
%                       matPi: matrix (banks x banks x seniorities) of
%                       relative liabilities
%               Outputs: conversion matrix (banks x banks x numK)
% * vecLambdaB: vector (banks x 1) of bail-in thresholds
% * vecLambdaR: vector (banks x 1) of recapitallization thresholds
% * numSimulations: number of simulation paths for Monte-Carlo
% * dT (optional): time increment. Default 1
%
% *Outputs*
%
% * matValuations: (banks x seniorities x numSimulations) array of payoff 
% matrices under different simulation paths. Payoffs are not discounted 
% yet!
%


function matValuations = calcValuationsBailIn(vecE,vecMu,matCovariance,matL,matTheta,numK,funConversion,vecLambdaB,vecLambdaR,numSimulations,dT)

%% Get inputs & Declarations
if nargin < 11 ; dT = 1; end

numBanks = length(vecE);
vecVariance = diag(matCovariance);
numSeniority = size(matL);
if length(numSeniority) == 3
    numSeniority = numSeniority(3);
else
    numSeniority = 1;
end

%% Simulate external asset returns
L = chol(matCovariance,'lower');
dW = L*randn(numBanks,numSimulations);
matE = repmat(vecE, 1, numSimulations);
matE = matE.*exp((vecMu-0.5*vecVariance)*dT+sqrt(dT)*dW);

%% Compute payoffs for each simulation path
matValuations = nan(numBanks,numSeniority,numSimulations);

for i = 1:numSimulations
    vecE_Final = matE(:,i);
    matPayments = calcPaymentsBailIn(vecE_Final,matL,matTheta,numK,funConversion,vecLambdaB,vecLambdaR);
    matValuations(:,:,i) = matPayments;
end

end