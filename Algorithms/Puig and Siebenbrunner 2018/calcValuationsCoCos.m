%% calcValuationsCoCos
% Monte-Carlo simulation of payoffs under geometric brownian motion for
% external assets. Uses the Puig/Siebenbrunner algorithm with contingent
% convertible debt to compute payoffs.
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
% * matConversion: conversion matrix (banks x banks), determines the
% amount of shares received in case of conversion
% * vecLambdaC: vector (banks x 1) of CoCo thresholds
% * vecF: vector (banks x 1) of conversion fractions, the fraction of debt to be converted
% if the CoCo is triggered
% * numSimulations: number of simulation paths for Monte-Carlo
% * dT (optional): time increment. Default 1
%
% *Outputs*
%
% * matValuations: (banks x seniorities x numSimulations) array of payoff 
% matrices under different simulation paths. Payoffs are not discounted 
% yet!
%


function matValuations = calcValuationsCoCos(vecE,vecMu,matCovariance,matL,matTheta,numK,matConversion,vecLambdaC,vecF,numSimulations,dT)

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
    matPayments = calcPaymentsCoCos(vecE_Final,matL,matTheta,numK,matConversion,vecLambdaC,vecF);
    matValuations(:,:,i) = matPayments;
end

end