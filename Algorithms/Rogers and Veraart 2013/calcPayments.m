%% calcPayments
% Computes a clearing payment vector for a given financial system, using
% the Eisenberg/Noe methodology. If liquidation haircuts are supplied, the
% Rogers/Veraart extension of the E/N algorithm is used
%
% Rogers, L. C., & Veraart, L. A. (2013). Failure and rescue in an
% interbank network. Management Science, 59(4), 882-898.
%
% *Inputs*
%
% * vecE: vector (banks x 1) of other assets
% * matL: matrix (banks x banks) of interbank claims
% * dblAlpha (optional): recovery value on other assets (default 1)
% * dblBeta (optional): recovery value on interbank assets (default 1)
% * strCalculationMethod (optional): 'Standard' for Eisenberg/Noe
% algorithm (default), 'Iterate' for iterative solution
%
% *Outputs*
%
% * vecP: clearing payment vector
% * vecEAfterContagion: value of other assets after haircuts
% * vecEquityAfterContagion: equity values after haircuts and contagion
% 
% Author: Christoph Siebenbrunner
% Last modified: 11.07.2015
%

function vecPayments = calcPayments(vecE,matL,dblAlpha,dblBeta,strCalculationMethod)

%% Get inputs & Declarations
if nargin < 5 ; strCalculationMethod = 'Standard'; end
if nargin < 4 ; dblBeta = 1; end
if nargin < 3 ; dblAlpha = 1; end

%%%
% Define E/N variables
numBanks = length(vecE);

vecPbar = sum(matL,2);
matPi = matL ./ repmat(vecPbar,1,numBanks);
matPi(isnan(matPi)) = 0;

%% Compute clearing payment vector
% Use Eisenberg/Noe-algorithm
if strcmp(strCalculationMethod,'Standard')
    vecPayments = vecPbar;
    blnLoop = true;

    mat = (eye(numBanks) - dblBeta*matPi');
    vecEquity = dblAlpha * vecE + dblBeta * matPi' * vecPbar  - vecPbar;
    while blnLoop
        posDefaulted=matPi'*vecPayments+vecE < vecPbar;
        vecPayments(posDefaulted)=mat(posDefaulted,posDefaulted)\vecEquity(posDefaulted)+vecPayments(posDefaulted);
        posDefaultedNew = matPi'*vecPayments+vecE < vecPbar;
        blnLoop = any(posDefaulted~=posDefaultedNew);
        if any(vecPayments<0)
            blnLoop = false;
            strCalculationMethod = 'Iterate';
            warning('calcPayments: No convergence in E/N-algorithm, iteration will be used. Check regularity conditions!')
        end
    end
end

%%%
% Use iterative algorithm (only if requested or E/N didn't converge)
if strcmp(strCalculationMethod,'Iterate')
    vecPayments = vecPbar;
    blnLoop = true;
    numMaxIterations = 500;
    iIteration = 1;
    
    while blnLoop
        vecPnew = vecPayments;
        posDefaulted = vecE + matPi' * vecPayments  < vecPbar; 
        vecLiquidatedValue = max(0,dblAlpha * vecE + dblBeta * matPi' * vecPayments);
        vecPnew(posDefaulted) = vecLiquidatedValue(posDefaulted);
        blnLoop = norm(abs(vecPnew-vecPayments))>0.0001;
        if iIteration == numMaxIterations
            warning('No convergence!')
            blnLoop = false;
        end
        vecPayments = vecPnew;
        iIteration = iIteration + 1;
    end
elseif ~strcmp(strCalculationMethod,'Standard')
    error('calcPayments: Unknown method - must be Standard or Iterate')
end

end