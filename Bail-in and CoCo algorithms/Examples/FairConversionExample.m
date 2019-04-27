% fair conversion example

%% Declarations
vecE = [0;50;10];
matL = [0 20 10;5 0 24;10 0 0];
matL(:,:,2) = [5 0 10;15 0 24;10 0 0];
matTheta = [0 0 0;0 0 0;0 0 0];
numK = 1;
vecBailIn = [15;0;0];
dblGamma = 0.5;

%% Auxiliary variables
matPi = matL;
matPbar = sum(matL,2);
for s=1:numSeniority
    matPi(:,:,s) = matL(:,:,s) ./ repmat(matPbar(:,s),1,numBanks);
end
matPi(isnan(matPi)) = 0;


%% Conversion 

funConversion = fairConversionFun(dblGamma);
matConversion = funConversion(vecBailIn, vecEquity, matPi);