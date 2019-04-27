%% Declarations
vecE = [0;50;10];
matL = [0 20 10;5 0 24;10 0 0];
matL(:,:,2) = [0 0 15;15 0 24;10 0 0];
matConversion = zeros(3); % write down 
matTheta = [0 0.1 0;0 0 0.15;0 0 0];
numK = 2;
vecLambdaC = zeros(3,1) + 0.05;
vecF = [0.5;0.5;0.5];



%% Conversion 
[matP, vecEquity, matTheta_New, matL_New] = calcElsingerCoCos(vecE,matL,matTheta,numK,matConversion,vecLambdaC,vecF);
