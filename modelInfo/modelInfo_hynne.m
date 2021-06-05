% Model/system info
modelInfo.S = sparse([[0   0   0   0   0   0   0   0   0   0   1  -1   0   0   0   0   0 -59   0   0   0   0   0   0;
                       0  59  -1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0;
                       0   0   0   0   0   0   0   0   0   0   0   0   0   0   1 -59   0   0   0   0   0   0   0   0;
                       0   0   0   0   0   0   0   0   0   0   0   1 -59   0   0   0   0   0   0   0   0   0   0   0;
                       0   0   1   0   1   0   0   0  -1  -1   0   0   0   0   0   0   0   0   0   0   0   1   1   2;
                       0   0   0   0   0   1   1  -1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0;
                       0   0   0   0   0   1  -1   0   0   0   0   0   0   0  -1   0   0   0   0   0   0   0   0   0;
                       0   0   1  -1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  -1   0   0;
                       0   0   0   0   0   0   0  -1   0   0   0   1   0   0   1   0   0   0   0   0   0   0   0   0;
                       0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1  -1  -1   0   0   0   0;
                       0   0   0   0   0   0   0   0   1  -1   0   0   0   0   0   0   0   0   0   0   0   0   0   0;
                       0   0   0   1  -1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0;
                       0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  -1   1   0   0   0;
                       0   0   0   0   0   0   0   0   0   1  -1   0   0   0   0   0   0   0   0   0   0   0   0   0;
                       0   0   0   0   0   0   0   0   0   0   0   0   1  -1   0   0   0   0   0   0   0   0   0   0;
                       0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1  -1   0   0   0   0   0   0   0;
                       1  -1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0;
                       0   0   0   0   0   0   0   1  -1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0;
                       0   0   0   0   1  -1   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0;
                       0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0  -1;
                       0   0  -1   0  -1   0   0   0   1   1   0   0   0   0   0   0   0   0   0   0   0  -1  -1  -1;
                       0   0   0   0   0   0   0   1   0   0   0  -1   0   0  -1   0   0   0   0   0   0   0   0   0]]);
                   
                     
                   
modelInfo.xBounds = [zeros(22,1) inf(22,1)];
                 
% Define flux bounds on the original model;  LB = -inf --> reversible
modelInfo.vBounds = [-inf(24,1) inf(24,1)];
% Set irreversible reactions
modelInfo.vBounds([1 3 5 10 11 12 14 15 17 19 20 22 23],1) = 0;

