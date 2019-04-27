%% calcPayments
% Computes a clearing payment vector for a given financial system, using
% the methodology of Eisenberg and Noe (2001).
%
% Eisenberg, L., & Noe, T. H. (2001). Systemic risk in financial systems. 
% Management Science, 47(2), 236-249.
%
% *Inputs*
%
% * vecE: vector (banks x 1) of other assets
% * matL: matrix (banks x banks) of interbank claims
% * strCalculationMethod (optional): 'Standard' for Eisenberg/Noe
% algorithm (default), 'Iterate' for iterative solution
% * numMaxIterations (optional): maximum number of iterations. Default 500
% * dblTolerance (optional): convergence threshold for iterative solution. 
% Default 0.0001
%
% *Outputs*
%
% * vecPayments: clearing payment vector
%

function vecPayments = calcPayments(vecE,matL,strCalculationMethod,numMaxIterations,dblTolerance)

%% Get inputs & Declarations
if nargin < 5 ; dblTolerance = 0.0001; end
if nargin < 4 ; numMaxIterations = 500; end
if nargin < 3 ; strCalculationMethod = 'Standard'; end

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

    mat = (eye(numBanks) - matPi');
    vecEquity = vecE + matPi' * vecPbar  - vecPbar;
    while blnLoop
        posDefaulted=matPi'*vecPayments+vecE < vecPbar;
        vecPayments(posDefaulted)=mat(posDefaulted,posDefaulted)\vecEquity(posDefaulted)+vecPbar(posDefaulted);
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
    iIteration = 1;
    
    while blnLoop
        vecPnew = vecPayments;
        posDefaulted = vecE + matPi' * vecPayments  < vecPbar; 
        vecLiquidatedValue = max(0,vecE + matPi' * vecPayments);
        vecPnew(posDefaulted) = vecLiquidatedValue(posDefaulted);
        blnLoop = norm(abs(vecPnew-vecPayments))>dblTolerance;
        if iIteration >= numMaxIterations
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