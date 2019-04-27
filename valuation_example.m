%% Declarations
vecE = [10;50;10]/10;

vecMu = [0.001;0.001;0.001];
matSigma = [0.1 0.005 -0.005;0.005 0.015 0.005;-0.005 0.005 0.1];
numTimeSteps = 100;
numSimulations = 1000;

matL = [0 20 10;5 0 24;10 0 0];
matL(:,:,2) = [0 0 15;15 0 24;10 0 0];
funConversion = @(x,y,z) zeros(size(x,1), size(x,1), size(x,2));

matTheta = [0 0.1 0;0 0 0.15;0 0 0];
numK = 1;
vecLambdaB = zeros(length(vecE),1) + 0.05;
vecLambdaR = zeros(length(vecE),1) + 0.2;

%% Create auxiliary variables
numSeniority = size(matL);
numSeniority = numSeniority(3);
numBanks = length(vecE);
matPi = matL;
matPbar = zeros(numBanks,numSeniority);
for s=1:numSeniority
    matPbar(:,s) = matPbar(:,s) + sum(matL(:,:,s),2);
end
for s=1:numSeniority
    matPi(:,:,s) = matL(:,:,s) ./ repmat(matPbar(:,s),1,numBanks);
end
matPi(isnan(matPi)) = 0;

claims = zeros(3,1);
liabilities = sum(sum(matL,3),2);
for i = 1:2
    claims = claims + matPi(:,:,i)'*matPbar(:,i);
end
equity = vecE + claims - liabilities;
total_assets = vecE + claims;
leverage = equity ./ total_assets;

%% Start simulations
matValuations = zeros(numBanks,numSimulations);
matClaims = zeros(numBanks,numSimulations);
vecE_Trajectories = GeomBrownianMotion_Multi(numSimulations,...
    1:numTimeSteps,1,vecMu,matSigma,vecE);

vecE_Final_Simulations = squeeze(vecE_Trajectories(:,:,end));

for i = 1:numSimulations
    vecE_Final = vecE_Final_Simulations(:,i);
    [matP, vecEquity, matTheta_New] = calcElsingerBailIn(vecE_Final,matL,matTheta,numK,funConversion,vecLambdaB,vecLambdaR);
    vecClaimsValue = zeros(numBanks,1);
    for j = 1:numSeniority
        vecClaimsValue = vecClaimsValue + matPi(:,:,j)'*matP(:,j);
    end
    vecParticipationValue = max(0,(matTheta_New - matTheta)) * max(0,vecEquity);
    matValuations(:,i) = vecClaimsValue + vecParticipationValue;
    matClaims(:,i) = vecClaimsValue;
end