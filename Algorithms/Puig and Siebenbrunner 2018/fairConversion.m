%% Fair Bail-in Conversion 
% Computes a fair conversion matrix.
%
% *Inputs*
%
%  matBailIn: matrix (banks x numK) of bail-in amount per bank and seniority class
%  dblGamma: parameter for fair conversion, as described in Puig and Siebenbrunner 2018. 
%  vecEquity: equity vector (banks x 1)
%  matPi: matrix (banks x banks x seniorities) ofrelative liabilities
%
%  *Outputs* 
%
%  conversion matrix (banks x banks x numK)
% 
% Author: Christoph Siebenbrunner and Matias Puig
% Last modified: 03.09.2018
%

function matConversion = fairConversion(matBailIn, dblGamma, vecEquity, matPi)
    numBanks = size(matBailIn,1);
    numK = size(matBailIn,2);
    matConversion = zeros(numBanks, numBanks, numK);
    numSeniority = size(matPi,3);
    for i=1:numBanks
        for j= 1:numBanks
            for k = 1:numK
                if vecEquity(j)> 0
                    matConversion(j,i,k) = matPi(j,i,numSeniority - numK + k)*matBailIn(j,k) / (vecEquity(j) + sum(matBailIn(j,:)));
                end
                if vecEquity(j) <= 0
                    if sum(matBailIn(j,:))>0
                        matConversion(j,i,k) = dblGamma * matPi(j,i,numSeniority - numK + k)*matBailIn(j,k) / sum(matBailIn(j,:));
                    end
                    if sum(matBailIn(j,:))==0
                        matConversion(j,i,k) = 0;
                    end
                end
            end
        end
    end
end


