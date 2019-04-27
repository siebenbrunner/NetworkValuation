%% calcElsingerOneLayer
% Computes a clearing payment vector for a given financial system, using
% the Elsinger (2009) methodology without seniority structure
%
% *Inputs*
%
% * vecE: vector (banks x 1) of other assets
% * matL: matrix (banks x banks) of interbank claims
% * matTheta: matrix (banks x banks) of relative interbank holdings
%
% *Outputs*
%
% * vecP: clearing payment vector
% * vecEquity: equity values after contagion
% 
% Author: Christoph Siebenbrunner and Matias Puig
% Last modified: 18.06.2018
%


function [vecP,vecEquity] = calcElsingerOneLayer(vecE,matL,matTheta)

%%%
% Variable initialisations
vecPbar = sum(matL,2);
matPi = matL ./ repmat(vecPbar,1,length(matL(:,1)));
matPi(isnan(matPi)) = 0;

% Convergence parameters
dblPrecision = max(vecPbar)/100000;
numMaxIterations=100;

%% Loop to find fixed point
% Initialise loop variables
vecP = vecPbar;
blnLoop = true;
numIterations=0;

while blnLoop
    vecP_old = vecP;
    
    %%%
    % Compute new equity value
    vecEquity = calcEquityValue(vecE,matPi,vecP,vecPbar,matTheta);
    posDefaulted = vecEquity < 0;
    
    %%%
    % Compute new p vector
    vecAux = matPi'*vecP+vecE+matTheta'*max(vecEquity,0);
    vecP(posDefaulted) = max(0,vecAux(posDefaulted));
    
    %%%
    % Check for convergence
    blnLoop = norm(abs(vecP-vecP_old)) > dblPrecision;
    if numIterations>numMaxIterations
        blnLoop=false;
        disp('No convergence in calcPaymentsOneLayer');
    end
    numIterations=numIterations+1;
end

end

%% Sub-function calcEquityValue
function vecEquity = calcEquityValue(vecE,matPi,vecP,vecPbar,matTheta)

%% Declarations
% Convergence parameters
dblPrecision = max(vecPbar)/100000;
numMaxIterations=100;
%%%
% Variable initialisations
vecEquity=matPi'*vecP+vecE-vecPbar;

%% Loop to find fixed point
% Initialise loop variables
blnLoop=true;
numIterations=0;
if sum(sum(matTheta))~=0
    while blnLoop
        vecEquity_old=vecEquity;
        %%%
        % Compute new equity vector
        vecEquity = max(vecEquity,0);
        vecEquity = vecE + matPi'*vecP - vecPbar + matTheta'*vecEquity;
        
        %%%
        % Check for convergence
        blnLoop = norm(abs(vecEquity-vecEquity_old)) > dblPrecision;
        if numIterations>numMaxIterations
            blnLoop=false;
            disp('No convergence in calcEquityValue');
        end
        numIterations=numIterations+1;
    end
end

end
