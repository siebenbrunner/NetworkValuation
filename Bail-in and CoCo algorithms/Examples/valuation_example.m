%% Declarations
vecE = [10;50;10];
vecMu = [0.001;0.001;0.001];
rho = 0.2; % correlations, could be generalized to be a matrix
dT = 1;
vecVariance= [0.1;0.1;0.1];
matStd = diag(sqrt(vecVariance));
numSimulations = 5;
numBanks = 3;
matCovariance = matStd * (eye(numBanks) + rho*(ones(numBanks)-eye(numBanks))) * matStd;

matL = [0 20 10;5 0 24;10 0 0];
matL(:,:,2) = [0 0 15;15 0 24;10 0 0];
funConversion = @(x,y,z) zeros(size(x,1), size(x,1), size(x,2)); % write down 

matTheta = [0 0.1 0;0 0 0.15;0 0 0];
numK = 1;
vecLambdaB = zeros(length(vecE),1) + 0.05;
vecLambdaR = zeros(length(vecE),1) + 0.2;


%% Simulation of asset paths
L = chol(matCovariance,'lower');
dW = L*randn(numBanks,numSimulations);
matE = repmat(vecE, 1, numSimulations);
matE = matE.*exp((vecMu-0.5*vecVariance)*dT+sqrt(dT)*dW);


%% Create auxiliary variables
numSeniority = size(matL,3);
matPi = matL;
matPbar = zeros(numBanks,numSeniority);
for s=1:numSeniority
    matPbar(:,s) = matPbar(:,s) + sum(matL(:,:,s),2);
end
for s=1:numSeniority
    matPi(:,:,s) = matL(:,:,s) ./ repmat(matPbar(:,s),1,numBanks);
end
matPi(isnan(matPi)) = 0;

%% Start simulations
matValuations = zeros(numBanks,numSimulations);

for i = 1:numSimulations
    vecE_Final = matE(:,i);
    [matP, vecEquity, matTheta_New] = calcElsingerBailIn(vecE_Final,matL,matTheta,numK,funConversion,vecLambdaB,vecLambdaR);
    matValuationSim = sum(matP,2) + max(matTheta_New - matTheta,0) * max(vecEquity,0);
    matValuations(:,i) = matValuations(:,i) + matValuationSim;
end

matValuations
matValuationAverage = sum(matValuations,2)/numSimulations;