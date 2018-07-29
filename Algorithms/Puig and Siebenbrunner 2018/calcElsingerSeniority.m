%% calcElsingerSeniority
% Computes a clearing payment matrix for a given financial system, using
% the Elsinger 2009 methodology with seniority structure. 
%
% *Inputs*
%
% * vecE: vector (banks x 1) of other assets
% * matL: matrix (banks x banks x seniorities) of interbank claims
% * matTheta: matrix (banks x banks) of interbank holdings

%
% *Outputs*
%
% * matP: clearing payment matrix (banks * seniorities)
% * vecEquityAfterContagion: equity values after contagion
% * matTheta: matrix (banks x banks x seniorities) of interbank holdings
% * vecDefaultedBanks: boolean vector (banks x 1), 1 if bank has defaulted

% Authors: Dissertation candidate and supervisor 
% Last modified: 18.06.2018
%

function [matP, vecEquity, matTheta, vecDefaultedBanks] = calcElsingerSeniority(vecE,matL,matTheta)

% Define Elsinger 2009 (seniority) variables

numSeniority = size(matL);
numSeniority = numSeniority(3);
numBanks = length(vecE);
vecDefaultedBanks = false(numBanks,1);
vecH = ones(numBanks,1) * numSeniority;
blnloop = true;
matLH = zeros(numBanks);
matPi = matL;
matPbar = zeros(numBanks,numSeniority);
for s=1:numSeniority
    matPbar(:,s) = matPbar(:,s) + sum(matL(:,:,s),2);
end
for s=1:numSeniority
    matPi(:,:,s) = matL(:,:,s) ./ repmat(matPbar(:,s),1,numBanks);
end
matPi(isnan(matPi)) = 0;

%% Compute clearing payment vector
% Use Elsinger 2009 algorithm with seniority structure

while blnloop
    % define vector vecEH
    vecEH = vecE;
    for i=1:numBanks
        for j=1:numBanks
            for s = 1:(vecH(j)-1)
                vecEH(i) = vecEH(i) + matPi(j,i,s)*matPbar(j,s);
            end
        end
    end
    
    for i=1:numBanks
        for s = 1:(vecH(i)-1)
            vecEH(i) = vecEH(i) - matPbar(i,s);
        end
    end
    
    % define matrix matLH
    for i=1:numBanks
        for j = 1:numBanks
            matLH(i,j) = matL(i,j,vecH(i));
        end
    end
    
            
    [vecP,vecEquity] = calcElsinger(vecEH,matLH, matTheta);
    blnDefault = vecEquity < 0 & vecH>1; 
    vecDefaultedBanks = vecDefaultedBanks | blnDefault;
    
    if (blnDefault == zeros(numBanks,1))
        blnloop = false;
    end
    

    % update variables 
    vecH = vecH - blnDefault;

matP = zeros(numBanks, numSeniority);

for i=1:numBanks
    for s = 1:numSeniority
        if(vecH(i) > s) 
            matP(i,s) = matPbar(i,s);
        end
        if(vecH(i) == s)
            matP(i,s) = vecP(i);
        end
        if(vecH(i) < s)
            matP(i,s) = 0;
        end
    end
end

% compute equity vector 

vecEquity = vecE; 
for s=1:numSeniority
    vecEquity = vecEquity + matPi(:,:,s)'*matP(:,s) - matPbar(:,s);
end

vecEquity = vecEquity + max(0,matTheta' * vecEquity);  




end

