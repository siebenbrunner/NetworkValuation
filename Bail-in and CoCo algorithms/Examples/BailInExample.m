%% Declarations
vecE = [0;50;10];
matL = [0 20 10;5 0 24;10 0 0];
matL(:,:,2) = [0 0 15;15 0 24;10 0 0];
funConversionWriteDown = @(x,y,z) zeros(size(x,1), size(x,1), size(x,2)); % write down
dblGamma = 0.5;
funConversionFair = fairConversionFun(dblGamma);
matTheta = [0 0 0;0 0 0;0 0 0];
numK = 1;
vecLambdaB = zeros(length(vecE),1) + 0.05;
vecLambdaR = zeros(length(vecE),1) + 0.2;



%% Bai-in
[matP_WriteDown, vecEquity_WriteDown, matTheta_WriteDown, matL_WriteDown] = calcElsingerBailIn(vecE,matL,matTheta,numK,funConversionWriteDown,vecLambdaB,vecLambdaR);
[matP_FairConversion, vecEquity_FairConversion, matTheta_FairConversion, matL_FairConversion] = calcElsingerBailIn(vecE,matL,matTheta,numK,funConversionFair,vecLambdaB,vecLambdaR);


%% recapitalization levels

vecLambdaWriteDown = vecEquity_WriteDown ./ (vecEquity_WriteDown + sum(sum(matL_WriteDown,3),2));
vecLambdaFairConversion = vecEquity_FairConversion ./ (vecEquity_FairConversion + sum(sum(matL_FairConversion,3),2));
