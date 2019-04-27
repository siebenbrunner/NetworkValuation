%% calcElsingerBailIn
% Computes a clearing payment vector for a given financial system, using
% the Elsinger 2009 methodology with seniority structure and bail-in. 
%
% *Inputs*
%
% * vecE: vector (banks x 1) of other assets
% * vecLambdaB: vector (banks x 1) of bail-in thresholds
% * vecLambdaR: vector (banks x 1) of recapitallization thresholds
% * matL: matrix (banks x banks x seniorities) of interbank claims
% * matTheta: matrix (banks x banks) of interbank holdings
% * numK: integer, determines the number of seniority classes corresponding to bail-in 
% * funConversion: function that returns the desired conversion factor.
%               Inputs: MatBailIn - matrix (banks x numK) of bail-in amount per bank and seniority class
%                       vecEquity: equity vector (banks x 1)
%                       matPi: matrix (banks x banks x seniorities) of
%                       relative liabilities
%               Outputs: conversion matrix (banks x banks x numK)
%
% *Outputs*
%
% * matP: clearing payment matrix (banks * seniorities)
% * vecEquityAfterContagion: equity values after contagion
% * matTheta: matrix (banks x banks) of interbank holdings
% * matL: liabilities matrix (banks x banks x seniorities)
% * vecDefaultedBanks: boolean vector (banks x 1), 1 if bank has defaulted
% * vecBailedInBanks: boolean vector (banks x 1), 1 if bank has been
%                     bailed-in
% Authors: Matias Puig and Christoph Siebenbrunner
% Last modified: 03.09.2018
%

function [matP, vecEquity, matTheta, matL, vecDefaultedBanks, vecBailedInBanks] = calcElsingerBailIn(vecE,matL, matTheta, numK, funConversion, vecLambdaB, vecLambdaR)


% Define Elsinger 2009 (seniority) variables

numSeniority = size(matL);
numSeniority = numSeniority(3);
numBanks = length(vecE);
vecDefaultedBanks = false(numBanks,1);
vecBailedInBanks = false(numBanks,1);
blnLoop = true;
numIterations = 0;
matF = zeros(numBanks, numBanks, numSeniority); % gained shares 
vecBailInAble = zeros(numBanks,1);
matPi = matL;
matPbar = zeros(numBanks,numSeniority); 
matBailIn = zeros(numBanks, numK);


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

%% Compute clearing payment vector
% Use Elsinger 2009 algorithm with seniority structure

while blnLoop
     [matP, vecEquity, matTheta, vecDefaultedBanks] = calcElsingerSeniority(vecE,matL,matTheta);
     
     matPbarOld = matPbar;
     % bail-in-able liabilities total
     
     for s=(numSeniority -numK+1):numSeniority
         vecBailInAble = vecBailInAble + matPbar(:,s);
     end
     vecLambda = vecEquity ./ (vecEquity + sum(matPbar,2)); % vector of capital ratios 
     vecBailIn = max(0,sum(matPbar,2) - (1 - vecLambdaR).*(vecEquity + sum(matPbar,2))); 
     vecBailIn(vecLambda > vecLambdaB) = 0; % bail in only if lambda < lambdaB
     vecDefaultedBanks = vecDefaultedBanks | (vecBailIn>0 & vecBailInAble==0);
     vecBailIn = min(vecBailInAble, vecBailIn);
     vecBailedInBanks = vecBailedInBanks | (vecBailIn>0);
     
     % bail-in process
     
     for s=(numSeniority-numK+1):numSeniority
         vecJunior = zeros(numBanks,1);
         for s2=(s+1):numSeniority
             vecJunior = vecJunior + matPbar(:,s2);
         end
         for i=1:numBanks
             if vecBailIn(i) >= vecJunior(i)
                 matPbar(i,s) = max(0,(matPbar(i,s)- vecBailIn(i) + vecJunior(i))); 
                 matBailIn(i,s + numK - numSeniority) = matPbarOld(i,s) - matPbar(i,s);
             end
             
         end
     end
     
     matConversion = funConversion(matBailIn, vecEquity, matPi);
     for i=1:numBanks
             for j=1:numBanks
                  for s=(numSeniority -numK+1):numSeniority
                    matF(i,j,s) =  matConversion(i,j,s + numK - numSeniority);
                  end
             end
     end
         
     vecFbar = sum(sum(matF, 3),2);
     
     
     
     % update liabilities array
     
     for s=1:numSeniority
         matL(:,:,s) = matPi(:,:,s).*matPbar(:,s);
     end
     
     % update holding matrix
     
     for i=1:numBanks
         for j=1:numBanks
             matTheta(j,i) = (1 - vecFbar(i))*matTheta(j,i) + sum(matF(j,i,:));
         end
     end
     
     blnLoop = norm(sum(abs(matPbarOld-matPbar),2)) > dblPrecision;
    if numIterations>numMaxIterations
        blnLoop=false;
        disp('No convergence in calcElsingerBailIn');
    end
    numIterations=numIterations+1;
end


end

