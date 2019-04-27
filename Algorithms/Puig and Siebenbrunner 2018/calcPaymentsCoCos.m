%% calcPaymentsCoCos
% Computes a clearing payment vector for a given financial system, using
% the Elsinger 2009 methodology with seniority structure and Contingent Convertibles (CoCo's). 
%
% *Inputs*
%
% * vecE: vector (banks x 1) of other assets
% * matL: matrix (banks x banks x seniorities) of interbank claims
% * matTheta: matrix (banks x banks) of interbank holdings
% * numK: integer, determines the number of the seniority class corresponding to CoCo's 
% * matConversion: conversion matrix (banks x banks), determines the
% amount of shares received in case of conversion
% * vecLambdaC: vector (banks x 1) of CoCo thresholds
% * vecF: vector (banks x 1) of conversion fractions, the fraction of debt to be converted
% if the CoCo is triggered
%
% *Outputs*
%
% * matP: clearing payment matrix (banks * seniorities)
% * vecEquityAfterContagion: equity values after contagion
% * matTheta: matrix (banks x banks x seniorities) of interbank holdings
% * matL: liabilities matrix (banks x banks x seniorities)
% * vecDefaultedBanks: boolean vector (banks x 1), 1 if bank has defaulted
% * vecTriggeredBanks: boolean vector (banks x 1), 1 if bank CoCo's have
% triggered
%    
% Authors: Matias Puig and Christoph Siebenbrunner
% Last modified: 03.09.2018
%

function [matP, vecEquity, matTheta, matL, vecDefaultedBanks, vecTriggeredBanks] = calcPaymentsCoCos(vecE,matL, matTheta, numK, matConversion, vecLambdaC, vecF)


% Define Elsinger 2009 (seniority) variables

numSeniority = size(matL,3);
numBanks = length(vecE);
vecDefaultedBanks = false(numBanks,1);
blnLoop = true;
numIterations = 0;
matPi = matL;
matPbar = zeros(numBanks,numSeniority); 



for s=1:numSeniority
    matPbar(:,s) = matPbar(:,s) + sum(matL(:,:,s),2);
end
for s=1:numSeniority
    matPi(:,:,s) = matL(:,:,s) ./ repmat(matPbar(:,s),1,numBanks);
end
matPi(isnan(matPi)) = 0;

% Convergence parameters
dblPrecision = max(matPbar)/100000;
numMaxIterations=100;


vecTriggeredBanks = [0;0;0];  % banks which have triggered CoCos

%% Compute clearing payment vector
% Use Elsinger 2009 algorithm with seniority structure

while blnLoop
     [matP, vecEquity, vecDefaultedBanks] = calcElsinger(vecE,matL,matTheta);
     
     matPbarOld = matPbar;
     vecLambda = vecEquity ./ (vecEquity + sum(matPbar,2)); % vector of capital ratios 
     vecTriggeredIteration = vecLambda < vecLambdaC; % triggered CoCos in this iteration
     vecTriggeredIteration = vecTriggeredIteration &  not(vecTriggeredBanks); % only trigger if not previously triggered

     
     % conversion process
     vecCbar = sum(matConversion,2); % amount of shares converted by each bank
     for i=1:numBanks
        if vecTriggeredIteration(i) 
            matPbar(i,numK) = (1-vecF(i))*matPbar(i,numK);
            for j=1:numBanks
                 matTheta(i,j) = (1 - vecCbar(i))*matTheta(i,j) + matConversion(i,j);
            end
        end
             
     end
   
    % update liabilities array
     
    for s=1:numSeniority
        matL(:,:,s) = matPi(:,:,s).*matPbar(:,s);
    end
     
    blnLoop = norm(sum(abs(matPbarOld-matPbar),2)) > dblPrecision;
    if numIterations>numMaxIterations
        blnLoop=false;
        disp('No convergence in calcElsingerBailIn');
    end
    numIterations=numIterations+1;
    vecTriggeredBanks = vecTriggeredBanks | vecTriggeredIteration;


end

