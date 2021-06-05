function [RF knownMet] = solveRF_KEapproach_hynne(relative_concMatrix,timeVec,knownKinetics)
% Kinetic equations approach, assuming some or all kinetic equations for
% reactions are known in the system.
% Parameters and kinetics were extracted from Full-scale model of
% glycolysis in Saccharomyces cerevisiae, Hynne and Sorensen 2001.

params(1) =  1014.96; % V2f
params(2) =  1.7; % K2Glc
params(3) =  1; % P2
params(4) =  1.2; % K2IG6P
params(5) =  7.2; % K2IIG6P
params(6) =  1014.96; % V2r
params(7) =  51.7547; % V3m
params(8) =  0.37; % K3DGlc
params(9) =  0.1; % K3ATP
params(10) =  0; % K3Glc
params(11) =  496.042; % V4f
params(12) =  0.8; % K4G6P
params(13) =  0.15; % K4F6P
params(14) =  496.042; % V4r
params(15) =  0.13; % K4eq
params(16) =  45.4327; % V5m
params(17) =  0.021; % K5
params(18) =  0.15; % kappa5
params(19) =  2207.82; % V6f
params(20) =  0.3; %K6FBP
params(21) =  2; % K6DHAP
params(22) =  0.081; % K6eq
params(23) =  5; % ratio6
params(24) =  4; % K6GAP
params(25) =  10; % K6IGAP
params(26) =  116.365; % V7f
params(27) =  1.23; % K7DHAP
params(28) =  1.27; % K7GAP
params(29) =  116.365; % V7r
params(30) =  0.055; % K7eq
params(31) =  833.858; % V8f
params(32) =  0.6; % K8GAP
params(33) =  0.1; % K8NAD
params(34) =  0.01; % K8BPG
params(35) =  0.06; % K8NADH
params(36) =  833.858; % V8r
params(37) =  0.0055; % K8eq
params(38) =  443866; % PEPsynth_kf
params(39) =  1528.62; % PEPsynth_kr
params(40) =  343.096; % V10m
params(41) =  0.2; % K10PEP
params(42) =  0.17; % K10ADP
params(43) =  53.1328; % V11m
params(44) =  0.3; % K11
params(45) =  89.8023; % V12m
params(46) =  0.1; % K12NADH
params(47) =  0.71; % K12ACA
params(48) =  16.72; % k13
params(49) =  81.4797; % V15m
params(50) =  25; % K15DHAP
params(51) =  0.034; % K15INADH
params(52) =  0.13; % K15INAD
params(53) =  0.13; % K15NADH
params(54) =  1.9; % k16
params(55) =  24.7; % k18
params(56) =  0.00283828; % k20
params(57) =  2.25932; % k22
params(58) =  3.2076; % k23
params(59) =  432.9; % Adenylate_kinase_kf
params(60) =  133.333; % Adenylate_kinase_kr

params(61) =  59; % Yvol
params(62) =  0.048; % k0

params(63) =  24; % GlcX0 (default 18.5)
params(64) =  5.6; % CNX0

for t = 2:size(relative_concMatrix,1)
    relative_Vpool(t-1,:) = (relative_concMatrix(t,:)-relative_concMatrix(t-1,:))./(timeVec(t)-timeVec(t-1));
end

knownMet = findKnownMet(knownKinetics);

minRF = max(relative_concMatrix(:,knownMet))/50;

lb = minRF;
ub = 1000*ones(length(knownMet),1);
x0 = (1000-1).*rand(1,length(knownMet)) + 1;

if ~isempty(knownMet)
    [RF,res] = lsqnonlin(@(RF) dxdt_functions(RF,relative_Vpool,relative_concMatrix,params,knownKinetics,knownMet),x0,lb,ub);
else
    RF = 0;
    knownMet = 0;
end

end

function func = dxdt_functions(RF,relative_Vpool,relative_concMatrix,params,knownKinetics,knownMet)

nT = size(relative_Vpool,1);
RFnum = zeros(size(relative_Vpool,2),1);
RFnum(knownMet) = 1:length(knownMet);

func = [];
for t = 2:nT
    
    % dxdt1
    if all(ismember([11 12 18],knownKinetics))
        func(end+1) = ...
            abs(((params(43)*(relative_concMatrix(t,14)/RF(RFnum(14))))/(params(44) + (relative_concMatrix(t,14)/RF(RFnum(14)))))...
            -1*((params(45)*(relative_concMatrix(t,1)/RF(RFnum(1)))*(relative_concMatrix(t,22)/RF(RFnum(22))))/((params(46) + (relative_concMatrix(t,22)/RF(RFnum(22))))*(params(47) + (relative_concMatrix(t,1)/RF(RFnum(1))))))...
            -59*((params(55)/params(61))*((relative_concMatrix(t,1)/RF(RFnum(1))) - (relative_concMatrix(t,10)/RF(RFnum(10)))))...
            -relative_Vpool(t,1)/RF(RFnum(1)));
    end
    
    % dxdt2
    if all(ismember([2 3],knownKinetics))
        func(end+1) = ...
            abs(59*((params(1)/params(61))*((relative_concMatrix(t,17)/RF(RFnum(17)))/params(2))/(1 + (relative_concMatrix(t,17)/RF(RFnum(17)))/params(2) + ((params(3)*(relative_concMatrix(t,17)/RF(RFnum(17)))/params(2) + 1)/(params(3)*(relative_concMatrix(t,2)/RF(RFnum(2)))/params(2) + 1)) * (1 + (relative_concMatrix(t,2)/RF(RFnum(2)))/params(2) + (relative_concMatrix(t,8)/RF(RFnum(8)))/params(4) + ((relative_concMatrix(t,2)/RF(RFnum(2)))*(relative_concMatrix(t,8)/RF(RFnum(8))))/(params(2)*params(5))))...
                - (params(6)/params(61))*((relative_concMatrix(t,2)/RF(RFnum(2)))/params(2))/(1 + (relative_concMatrix(t,2)/RF(RFnum(2)))/params(2) + ((params(3)*(relative_concMatrix(t,2)/RF(RFnum(2)))/params(2) + 1)/(params(3)*(relative_concMatrix(t,17)/RF(RFnum(17)))/params(2) + 1)) * (1 + (relative_concMatrix(t,2)/RF(RFnum(2)))/params(2)) + (relative_concMatrix(t,8)/RF(RFnum(8)))/params(4) + ((relative_concMatrix(t,2)/RF(RFnum(2)))*(relative_concMatrix(t,8)/RF(RFnum(8))))/(params(2)*params(5))))...
            -1*((params(7)*(relative_concMatrix(t,21)/RF(RFnum(21)))*(relative_concMatrix(t,2)/RF(RFnum(2))))/(params(8)*params(9) + params(10)*(relative_concMatrix(t,21)/RF(RFnum(21))) + params(9)*(relative_concMatrix(t,2)/RF(RFnum(2))) + (relative_concMatrix(t,2)/RF(RFnum(2)))*(relative_concMatrix(t,21)/RF(RFnum(21)))))...
            -relative_Vpool(t,2)/RF(RFnum(2)));
    end
    
    % dxdt3
    if all(ismember([15 16],knownKinetics))
        func(end+1) = ...
            abs(((params(49)*(relative_concMatrix(t,7)/RF(RFnum(7))))/(params(50)*(1 + (params(51)/(relative_concMatrix(t,22)/RF(RFnum(22))))*(1 + (relative_concMatrix(t,9)/RF(RFnum(9)))/params(52))) + (relative_concMatrix(t,7)/RF(RFnum(7)))*(1 + (params(53)/(relative_concMatrix(t,22)/RF(RFnum(22))))*(1 + (relative_concMatrix(t,9)/RF(RFnum(9)))/params(52)))))...
            -59*((params(54)/params(61))*((relative_concMatrix(t,3)/RF(RFnum(3))) - (relative_concMatrix(t,16)/RF(RFnum(16)))))...
            -relative_Vpool(t,3)/RF(RFnum(3)));
    end
    
    % dxdt4    
    if all(ismember([12 13],knownKinetics))
        func(end+1) = ...
            abs(((params(45)*(relative_concMatrix(t,1)/RF(RFnum(1)))*(relative_concMatrix(t,22)/RF(RFnum(22))))/((params(46) + (relative_concMatrix(t,22)/RF(RFnum(22))))*(params(47) + (relative_concMatrix(t,1)/RF(RFnum(1))))))...
            -59*((params(48)/params(61))*((relative_concMatrix(t,4)/RF(RFnum(4))) - (relative_concMatrix(t,15)/RF(RFnum(15)))))...
            -relative_Vpool(t,4)/RF(RFnum(4)));
    end
    
    % dxdt5
    if all(ismember([3 5 9 10 22 23 24],knownKinetics))
        func(end+1) = ...
            abs(((params(7)*(relative_concMatrix(t,21)/RF(RFnum(21)))*(relative_concMatrix(t,2)/RF(RFnum(2))))/(params(8)*params(9) + params(10)*(relative_concMatrix(t,21)/RF(RFnum(21))) + params(9)*(relative_concMatrix(t,2)/RF(RFnum(2))) + (relative_concMatrix(t,2)/RF(RFnum(2)))*(relative_concMatrix(t,21)/RF(RFnum(21)))))...
            +1*((params(16)*(relative_concMatrix(t,12)/RF(RFnum(12)))^2)/(params(17)*(1 + params(18)*((relative_concMatrix(t,21)/RF(RFnum(21)))/(relative_concMatrix(t,20)/RF(RFnum(20))))*((relative_concMatrix(t,21)/RF(RFnum(21)))/(relative_concMatrix(t,20)/RF(RFnum(20))))) + (relative_concMatrix(t,12)/RF(RFnum(12)))^2))...
            -1*(params(38)*(relative_concMatrix(t,18)/RF(RFnum(18)))*(relative_concMatrix(t,5)/RF(RFnum(5))) - params(39)*(relative_concMatrix(t,11)/RF(RFnum(11)))*(relative_concMatrix(t,21)/RF(RFnum(21))))...
            -1*((params(40)*(relative_concMatrix(t,5)/RF(RFnum(5)))*(relative_concMatrix(t,11)/RF(RFnum(11))))/((params(41) + (relative_concMatrix(t,11)/RF(RFnum(11))))*(params(42) + (relative_concMatrix(t,5)/RF(RFnum(5))))))...
            +1*(params(57)*(relative_concMatrix(t,21)/RF(RFnum(21)))*(relative_concMatrix(t,8)/RF(RFnum(8))))...
            +1*(params(58)*(relative_concMatrix(t,21)/RF(RFnum(21))))...
            +2*(params(59)*(relative_concMatrix(t,21)/RF(RFnum(21)))*(relative_concMatrix(t,20)/RF(RFnum(20))) - params(60)*(relative_concMatrix(t,5)/RF(RFnum(5)))*(relative_concMatrix(t,5)/RF(RFnum(5))))...
            -relative_Vpool(t,5)/RF(RFnum(5)));
    end
    
    % dxdt6
    if all(ismember([6 7 8],knownKinetics))
        func(end+1) = ...
            abs(((params(19)*(relative_concMatrix(t,19)/RF(RFnum(19))))/(params(20) + (relative_concMatrix(t,19)/RF(RFnum(19))) + ((relative_concMatrix(t,6)/RF(RFnum(6)))*params(21)*params(19))/(params(22)*params(19)*params(23)) + ((relative_concMatrix(t,7)/RF(RFnum(7)))*params(24)*params(19))/(params(22)*params(19)*params(23)) + (relative_concMatrix(t,19)/RF(RFnum(19)))*(relative_concMatrix(t,6)/RF(RFnum(6)))/params(25) + ((relative_concMatrix(t,6)/RF(RFnum(6)))*(relative_concMatrix(t,7)/RF(RFnum(7)))*params(19))/(params(22)*params(19)*params(23)))...
                - ((params(19)*(relative_concMatrix(t,6)/RF(RFnum(6)))*(relative_concMatrix(t,7)/RF(RFnum(7))))/params(22))/(params(20) + (relative_concMatrix(t,19)/RF(RFnum(19))) + ((relative_concMatrix(t,6)/RF(RFnum(6)))*params(21)*params(19))/(params(22)*params(19)*params(23)) + ((relative_concMatrix(t,7)/RF(RFnum(7)))*params(24)*params(19))/(params(22)*params(19)*params(23)) + (relative_concMatrix(t,19)/RF(RFnum(19)))*(relative_concMatrix(t,6)/RF(RFnum(6)))/params(25) + ((relative_concMatrix(t,6)/RF(RFnum(6)))*(relative_concMatrix(t,7)/RF(RFnum(7)))*params(19))/(params(22)*params(19)*params(23))))...
            +1*((params(26)*(relative_concMatrix(t,7)/RF(RFnum(7))))/(params(27) + (relative_concMatrix(t,7)/RF(RFnum(7))) + (params(27)/params(28))*(relative_concMatrix(t,6)/RF(RFnum(6)))) - (params(29)*(relative_concMatrix(t,6)/RF(RFnum(6)))/params(30))/(params(27) + (relative_concMatrix(t,7)/RF(RFnum(7))) + (params(27)/params(28))*(relative_concMatrix(t,6)/RF(RFnum(6)))))...
            -1*(((params(31)*(relative_concMatrix(t,6)/RF(RFnum(6)))*(relative_concMatrix(t,9)/RF(RFnum(9))))/params(32)/params(33))/((1 + (relative_concMatrix(t,6)/RF(RFnum(6)))/params(32) + (relative_concMatrix(t,18)/RF(RFnum(18)))/params(34))*(1 + (relative_concMatrix(t,9)/RF(RFnum(9)))/params(33) + (relative_concMatrix(t,22)/RF(RFnum(22)))/params(35)))...
                - ((params(36)*(relative_concMatrix(t,18)/RF(RFnum(18)))*(relative_concMatrix(t,22)/RF(RFnum(22))))/params(37)/params(32)/params(33))/((1 + (relative_concMatrix(t,6)/RF(RFnum(6)))/params(32) + (relative_concMatrix(t,18)/RF(RFnum(18)))/params(34))*(1 + (relative_concMatrix(t,9)/RF(RFnum(9)))/params(33) + (relative_concMatrix(t,22)/RF(RFnum(22)))/params(35))))...
            -relative_Vpool(t,6)/RF(RFnum(6)));
    end
    
    % dxdt7
    if all(ismember([6 7 15],knownKinetics))
        func(end+1) = ...
            abs(((params(19)*(relative_concMatrix(t,19)/RF(RFnum(19))))/(params(20) + (relative_concMatrix(t,19)/RF(RFnum(19))) + ((relative_concMatrix(t,6)/RF(RFnum(6)))*params(21)*params(19))/(params(22)*params(19)*params(23)) + ((relative_concMatrix(t,7)/RF(RFnum(7)))*params(24)*params(19))/(params(22)*params(19)*params(23)) + (relative_concMatrix(t,19)/RF(RFnum(19)))*(relative_concMatrix(t,6)/RF(RFnum(6)))/params(25) + ((relative_concMatrix(t,6)/RF(RFnum(6)))*(relative_concMatrix(t,7)/RF(RFnum(7)))*params(19))/(params(22)*params(19)*params(23)))...
                - ((params(19)*(relative_concMatrix(t,6)/RF(RFnum(6)))*(relative_concMatrix(t,7)/RF(RFnum(7))))/params(22))/(params(20) + (relative_concMatrix(t,19)/RF(RFnum(19))) + ((relative_concMatrix(t,6)/RF(RFnum(6)))*params(21)*params(19))/(params(22)*params(19)*params(23)) + ((relative_concMatrix(t,7)/RF(RFnum(7)))*params(24)*params(19))/(params(22)*params(19)*params(23)) + (relative_concMatrix(t,19)/RF(RFnum(19)))*(relative_concMatrix(t,6)/RF(RFnum(6)))/params(25) + ((relative_concMatrix(t,6)/RF(RFnum(6)))*(relative_concMatrix(t,7)/RF(RFnum(7)))*params(19))/(params(22)*params(19)*params(23))))...
            -1*((params(26)*(relative_concMatrix(t,7)/RF(RFnum(7))))/(params(27) + (relative_concMatrix(t,7)/RF(RFnum(7))) + (params(27)/params(28))*(relative_concMatrix(t,6)/RF(RFnum(6)))) - (params(29)*(relative_concMatrix(t,6)/RF(RFnum(6)))/params(30))/(params(27) + (relative_concMatrix(t,7)/RF(RFnum(7))) + (params(27)/params(28))*(relative_concMatrix(t,6)/RF(RFnum(6)))))...
            -1*((params(49)*(relative_concMatrix(t,7)/RF(RFnum(7))))/(params(50)*(1 + (params(51)/(relative_concMatrix(t,22)/RF(RFnum(22))))*(1 + (relative_concMatrix(t,9)/RF(RFnum(9)))/params(52))) + (relative_concMatrix(t,7)/RF(RFnum(7)))*(1 + (params(53)/(relative_concMatrix(t,22)/RF(RFnum(22))))*(1 + (relative_concMatrix(t,9)/RF(RFnum(9)))/params(52)))))...
            -relative_Vpool(t,7)/RF(RFnum(7)));
    end
    
    % dxdt8
    if all(ismember([3 4 22],knownKinetics))
        func(end+1) = ...
            abs(((params(7)*(relative_concMatrix(t,21)/RF(RFnum(21)))*(relative_concMatrix(t,2)/RF(RFnum(2))))/(params(8)*params(9) + params(10)*(relative_concMatrix(t,21)/RF(RFnum(21))) + params(9)*(relative_concMatrix(t,2)/RF(RFnum(2))) + (relative_concMatrix(t,2)/RF(RFnum(2)))*(relative_concMatrix(t,21)/RF(RFnum(21)))))...
            -1*((params(11)*(relative_concMatrix(t,8)/RF(RFnum(8))))/(params(12) + (relative_concMatrix(t,8)/RF(RFnum(8))) + (params(12)/params(13))*(relative_concMatrix(t,12)/RF(RFnum(12))))...
                - (params(14)*(relative_concMatrix(t,12)/RF(RFnum(12)))/params(15))/(params(12) + (relative_concMatrix(t,8)/RF(RFnum(8))) + (params(12)/params(13))*(relative_concMatrix(t,12)/RF(RFnum(12)))))...
            -1*(params(57)*(relative_concMatrix(t,21)/RF(RFnum(21)))*(relative_concMatrix(t,8)/RF(RFnum(8))))...
            -relative_Vpool(t,8)/RF(RFnum(8)));
    end
    
    % dxdt9
    if all(ismember([8 12 15],knownKinetics))
        func(end+1) = ...
            abs(-1*(((params(31)*(relative_concMatrix(t,6)/RF(RFnum(6)))*(relative_concMatrix(t,9)/RF(RFnum(9))))/params(32)/params(33))/((1 + (relative_concMatrix(t,6)/RF(RFnum(6)))/params(32) + (relative_concMatrix(t,18)/RF(RFnum(18)))/params(34))*(1 + (relative_concMatrix(t,9)/RF(RFnum(9)))/params(33) + (relative_concMatrix(t,22)/RF(RFnum(22)))/params(35)))...
                - ((params(36)*(relative_concMatrix(t,18)/RF(RFnum(18)))*(relative_concMatrix(t,22)/RF(RFnum(22))))/params(37)/params(32)/params(33))/((1 + (relative_concMatrix(t,6)/RF(RFnum(6)))/params(32) + (relative_concMatrix(t,18)/RF(RFnum(18)))/params(34))*(1 + (relative_concMatrix(t,9)/RF(RFnum(9)))/params(33) + (relative_concMatrix(t,22)/RF(RFnum(22)))/params(35))))...
            +1*((params(45)*(relative_concMatrix(t,1)/RF(RFnum(1)))*(relative_concMatrix(t,22)/RF(RFnum(22))))/((params(46) + (relative_concMatrix(t,22)/RF(RFnum(22))))*(params(47) + (relative_concMatrix(t,1)/RF(RFnum(1))))))...
            +1*((params(49)*(relative_concMatrix(t,7)/RF(RFnum(7))))/(params(50)*(1 + (params(51)/(relative_concMatrix(t,22)/RF(RFnum(22))))*(1 + (relative_concMatrix(t,9)/RF(RFnum(9)))/params(52))) + (relative_concMatrix(t,7)/RF(RFnum(7)))*(1 + (params(53)/(relative_concMatrix(t,22)/RF(RFnum(22))))*(1 + (relative_concMatrix(t,9)/RF(RFnum(9)))/params(52)))))...
            -relative_Vpool(t,9)/RF(RFnum(9)));
    end
    
    % dxdt10
    if all(ismember([8 20],knownKinetics)) && sum(ismember([14 17 19],knownKinetics)) > 0
        func(end+1) = ...
            abs(((params(55)/params(61))*((relative_concMatrix(t,1)/RF(RFnum(1))) - (relative_concMatrix(t,10)/RF(RFnum(10)))))...
            -1*(params(62)*(relative_concMatrix(t,10)/RF(RFnum(10))))...
            -1*(params(56)*(relative_concMatrix(t,10)/RF(RFnum(10)))*(relative_concMatrix(t,13)/RF(RFnum(13))))...
            -relative_Vpool(t,10)/RF(RFnum(10)));
    end
    
    % dxdt11
    if all(ismember([9 10],knownKinetics))
        func(end+1) = ...
            abs((params(38)*(relative_concMatrix(t,18)/RF(RFnum(18)))*(relative_concMatrix(t,5)/RF(RFnum(5))) - params(39)*(relative_concMatrix(t,11)/RF(RFnum(11)))*(relative_concMatrix(t,21)/RF(RFnum(21))))...
            -1*((params(40)*(relative_concMatrix(t,5)/RF(RFnum(5)))*(relative_concMatrix(t,11)/RF(RFnum(11))))/((params(41) + (relative_concMatrix(t,11)/RF(RFnum(11))))*(params(42) + (relative_concMatrix(t,5)/RF(RFnum(5))))))...
            -relative_Vpool(t,11)/RF(RFnum(11)));
    end
    
    % dxdt12
    if all(ismember([4 5],knownKinetics))
        func(end+1) = ...
            abs(((params(11)*(relative_concMatrix(t,8)/RF(RFnum(8))))/(params(12) + (relative_concMatrix(t,8)/RF(RFnum(8))) + (params(12)/params(13))*(relative_concMatrix(t,12)/RF(RFnum(12))))...
                - (params(14)*(relative_concMatrix(t,12)/RF(RFnum(12)))/params(15))/(params(12) + (relative_concMatrix(t,8)/RF(RFnum(8))) + (params(12)/params(13))*(relative_concMatrix(t,12)/RF(RFnum(12)))))...
            -1*((params(16)*(relative_concMatrix(t,12)/RF(RFnum(12)))^2)/(params(17)*(1 + params(18)*((relative_concMatrix(t,21)/RF(RFnum(21)))/(relative_concMatrix(t,20)/RF(RFnum(20))))*((relative_concMatrix(t,21)/RF(RFnum(21)))/(relative_concMatrix(t,20)/RF(RFnum(20))))) + (relative_concMatrix(t,12)/RF(RFnum(12)))^2))...
            -relative_Vpool(t,12)/RF(RFnum(12)));
    end
    
    % dxdt13
    if all(ismember([20 21],knownKinetics))
        func(end+1) = ...
            abs(-1*(params(56)*(relative_concMatrix(t,10)/RF(RFnum(10)))*(relative_concMatrix(t,13)/RF(RFnum(13))))...
            +1*(params(62)*(params(64) - (relative_concMatrix(t,13)/RF(RFnum(13)))))...
            -relative_Vpool(t,13)/RF(RFnum(13)));
    end
    
    % dxdt14
    if all(ismember([10 11],knownKinetics))
        func(end+1) = ...
            abs(((params(40)*(relative_concMatrix(t,5)/RF(RFnum(5)))*(relative_concMatrix(t,11)/RF(RFnum(11))))/((params(41) + (relative_concMatrix(t,11)/RF(RFnum(11))))*(params(42) + (relative_concMatrix(t,5)/RF(RFnum(5))))))...
            -1*((params(43)*(relative_concMatrix(t,14)/RF(RFnum(14))))/(params(44) + (relative_concMatrix(t,14)/RF(RFnum(14)))))...
            -relative_Vpool(t,14)/RF(RFnum(14)));
    end
    
    % dxdt15
    if all(ismember([13],knownKinetics)) && sum(ismember([14 17 19],knownKinetics)) > 0
        func(end+1) = ...
            abs(((params(48)/params(61))*((relative_concMatrix(t,4)/RF(RFnum(4))) - (relative_concMatrix(t,15)/RF(RFnum(15)))))...
            -1*(params(62)*(relative_concMatrix(t,15)/RF(RFnum(15))))...
            -relative_Vpool(t,15)/RF(RFnum(15)));
    end
    
    % dxdt16
    if all(ismember([16],knownKinetics)) && sum(ismember([14 17 19],knownKinetics)) > 0
        func(end+1) = ...
            abs(((params(54)/params(61))*((relative_concMatrix(t,3)/RF(RFnum(3))) - (relative_concMatrix(t,16)/RF(RFnum(16)))))...
            -1*(params(62)*(relative_concMatrix(t,16)/RF(RFnum(16))))...
            -relative_Vpool(t,16)/RF(RFnum(16)));
    end
    
    % dxdt17
    if all(ismember([1 2],knownKinetics))
        func(end+1) = ...
            abs((params(62)*(params(63)-(relative_concMatrix(t,17)/RF(RFnum(17)))))...
            -1*((params(1)/params(61))*((relative_concMatrix(t,17)/RF(RFnum(17)))/params(2))/(1 + (relative_concMatrix(t,17)/RF(RFnum(17)))/params(2) + ((params(3)*(relative_concMatrix(t,17)/RF(RFnum(17)))/params(2) + 1)/(params(3)*(relative_concMatrix(t,2)/RF(RFnum(2)))/params(2) + 1)) * (1 + (relative_concMatrix(t,2)/RF(RFnum(2)))/params(2) + (relative_concMatrix(t,8)/RF(RFnum(8)))/params(4) + ((relative_concMatrix(t,2)/RF(RFnum(2)))*(relative_concMatrix(t,8)/RF(RFnum(8))))/(params(2)*params(5))))...
                - (params(6)/params(61))*((relative_concMatrix(t,2)/RF(RFnum(2)))/params(2))/(1 + (relative_concMatrix(t,2)/RF(RFnum(2)))/params(2) + ((params(3)*(relative_concMatrix(t,2)/RF(RFnum(2)))/params(2) + 1)/(params(3)*(relative_concMatrix(t,17)/RF(RFnum(17)))/params(2) + 1)) * (1 + (relative_concMatrix(t,2)/RF(RFnum(2)))/params(2)) + (relative_concMatrix(t,8)/RF(RFnum(8)))/params(4) + ((relative_concMatrix(t,2)/RF(RFnum(2)))*(relative_concMatrix(t,8)/RF(RFnum(8))))/(params(2)*params(5))))...
            -relative_Vpool(t,17)/RF(RFnum(17)));
    end
    
    % dxdt18
    if all(ismember([8 9],knownKinetics))
        func(end+1) = ...
            abs((((params(31)*(relative_concMatrix(t,6)/RF(RFnum(6)))*(relative_concMatrix(t,9)/RF(RFnum(9))))/params(32)/params(33))/((1 + (relative_concMatrix(t,6)/RF(RFnum(6)))/params(32) + (relative_concMatrix(t,18)/RF(RFnum(18)))/params(34))*(1 + (relative_concMatrix(t,9)/RF(RFnum(9)))/params(33) + (relative_concMatrix(t,22)/RF(RFnum(22)))/params(35)))...
                - ((params(36)*(relative_concMatrix(t,18)/RF(RFnum(18)))*(relative_concMatrix(t,22)/RF(RFnum(22))))/params(37)/params(32)/params(33))/((1 + (relative_concMatrix(t,6)/RF(RFnum(6)))/params(32) + (relative_concMatrix(t,18)/RF(RFnum(18)))/params(34))*(1 + (relative_concMatrix(t,9)/RF(RFnum(9)))/params(33) + (relative_concMatrix(t,22)/RF(RFnum(22)))/params(35))))...
            -1*(params(38)*(relative_concMatrix(t,18)/RF(RFnum(18)))*(relative_concMatrix(t,5)/RF(RFnum(5))) - params(39)*(relative_concMatrix(t,11)/RF(RFnum(11)))*(relative_concMatrix(t,21)/RF(RFnum(21))))...
            -relative_Vpool(t,18)/RF(RFnum(18)));
    end
    
    % dxdt19
    if all(ismember([5 6],knownKinetics))
        func(end+1) = ...
            abs(((params(16)*(relative_concMatrix(t,12)/RF(RFnum(12)))^2)/(params(17)*(1 + params(18)*((relative_concMatrix(t,21)/RF(RFnum(21)))/(relative_concMatrix(t,20)/RF(RFnum(20))))*((relative_concMatrix(t,21)/RF(RFnum(21)))/(relative_concMatrix(t,20)/RF(RFnum(20))))) + (relative_concMatrix(t,12)/RF(RFnum(12)))^2))...
            -1*((params(19)*(relative_concMatrix(t,19)/RF(RFnum(19))))/(params(20) + (relative_concMatrix(t,19)/RF(RFnum(19))) + ((relative_concMatrix(t,6)/RF(RFnum(6)))*params(21)*params(19))/(params(22)*params(19)*params(23)) + ((relative_concMatrix(t,7)/RF(RFnum(7)))*params(24)*params(19))/(params(22)*params(19)*params(23)) + (relative_concMatrix(t,19)/RF(RFnum(19)))*(relative_concMatrix(t,6)/RF(RFnum(6)))/params(25) + ((relative_concMatrix(t,6)/RF(RFnum(6)))*(relative_concMatrix(t,7)/RF(RFnum(7)))*params(19))/(params(22)*params(19)*params(23)))...
                - ((params(19)*(relative_concMatrix(t,6)/RF(RFnum(6)))*(relative_concMatrix(t,7)/RF(RFnum(7))))/params(22))/(params(20) + (relative_concMatrix(t,19)/RF(RFnum(19))) + ((relative_concMatrix(t,6)/RF(RFnum(6)))*params(21)*params(19))/(params(22)*params(19)*params(23)) + ((relative_concMatrix(t,7)/RF(RFnum(7)))*params(24)*params(19))/(params(22)*params(19)*params(23)) + (relative_concMatrix(t,19)/RF(RFnum(19)))*(relative_concMatrix(t,6)/RF(RFnum(6)))/params(25) + ((relative_concMatrix(t,6)/RF(RFnum(6)))*(relative_concMatrix(t,7)/RF(RFnum(7)))*params(19))/(params(22)*params(19)*params(23))))...
            -relative_Vpool(t,19)/RF(RFnum(19)));
    end
    
    % dxdt20
    if all(ismember([24],knownKinetics))
        func(end+1) = ...
            abs(-1*(params(56)*(relative_concMatrix(t,10)/RF(RFnum(10)))*(relative_concMatrix(t,13)/RF(RFnum(13))))...
            -relative_Vpool(t,20)/RF(RFnum(20)));
    end
    
    % dxdt21
    if all(ismember([3 5 9 10 22 23 24],knownKinetics))
        func(end+1) = ...
            abs(-1*((params(7)*(relative_concMatrix(t,21)/RF(RFnum(21)))*(relative_concMatrix(t,2)/RF(RFnum(2))))/(params(8)*params(9) + params(10)*(relative_concMatrix(t,21)/RF(RFnum(21))) + params(9)*(relative_concMatrix(t,2)/RF(RFnum(2))) + (relative_concMatrix(t,2)/RF(RFnum(2)))*(relative_concMatrix(t,21)/RF(RFnum(21)))))...
            -1*((params(16)*(relative_concMatrix(t,12)/RF(RFnum(12)))^2)/(params(17)*(1 + params(18)*((relative_concMatrix(t,21)/RF(RFnum(21)))/(relative_concMatrix(t,20)/RF(RFnum(20))))*((relative_concMatrix(t,21)/RF(RFnum(21)))/(relative_concMatrix(t,20)/RF(RFnum(20))))) + (relative_concMatrix(t,12)/RF(RFnum(12)))^2))...
            +1*(params(38)*(relative_concMatrix(t,18)/RF(RFnum(18)))*(relative_concMatrix(t,5)/RF(RFnum(5))) - params(39)*(relative_concMatrix(t,11)/RF(RFnum(11)))*(relative_concMatrix(t,21)/RF(RFnum(21))))...
            +1*((params(40)*(relative_concMatrix(t,5)/RF(RFnum(5)))*(relative_concMatrix(t,11)/RF(RFnum(11))))/((params(41) + (relative_concMatrix(t,11)/RF(RFnum(11))))*(params(42) + (relative_concMatrix(t,5)/RF(RFnum(5))))))...
            -1*(params(57)*(relative_concMatrix(t,21)/RF(RFnum(21)))*(relative_concMatrix(t,8)/RF(RFnum(8))))...
            -1*(params(58)*(relative_concMatrix(t,21)/RF(RFnum(21))))...
            -1*(params(59)*(relative_concMatrix(t,21)/RF(RFnum(21)))*(relative_concMatrix(t,20)/RF(RFnum(20))) - params(60)*(relative_concMatrix(t,5)/RF(RFnum(5)))*(relative_concMatrix(t,5)/RF(RFnum(5))))...
            -relative_Vpool(t,21)/RF(RFnum(21)));
    end
    
    % dxdt22
    if all(ismember([8 12 15],knownKinetics))
        func(end+1) = ...
            abs((((params(31)*(relative_concMatrix(t,6)/RF(RFnum(6)))*(relative_concMatrix(t,9)/RF(RFnum(9))))/params(32)/params(33))/((1 + (relative_concMatrix(t,6)/RF(RFnum(6)))/params(32) + (relative_concMatrix(t,18)/RF(RFnum(18)))/params(34))*(1 + (relative_concMatrix(t,9)/RF(RFnum(9)))/params(33) + (relative_concMatrix(t,22)/RF(RFnum(22)))/params(35)))...
                - ((params(36)*(relative_concMatrix(t,18)/RF(RFnum(18)))*(relative_concMatrix(t,22)/RF(RFnum(22))))/params(37)/params(32)/params(33))/((1 + (relative_concMatrix(t,6)/RF(RFnum(6)))/params(32) + (relative_concMatrix(t,18)/RF(RFnum(18)))/params(34))*(1 + (relative_concMatrix(t,9)/RF(RFnum(9)))/params(33) + (relative_concMatrix(t,22)/RF(RFnum(22)))/params(35))))...
            -1*((params(45)*(relative_concMatrix(t,1)/RF(RFnum(1)))*(relative_concMatrix(t,22)/RF(RFnum(22))))/((params(46) + (relative_concMatrix(t,22)/RF(RFnum(22))))*(params(47) + (relative_concMatrix(t,1)/RF(RFnum(1))))))...
            -1*((params(49)*(relative_concMatrix(t,7)/RF(RFnum(7))))/(params(50)*(1 + (params(51)/(relative_concMatrix(t,22)/RF(RFnum(22))))*(1 + (relative_concMatrix(t,9)/RF(RFnum(9)))/params(52))) + (relative_concMatrix(t,7)/RF(RFnum(7)))*(1 + (params(53)/(relative_concMatrix(t,22)/RF(RFnum(22))))*(1 + (relative_concMatrix(t,9)/RF(RFnum(9)))/params(52)))))...
            -relative_Vpool(t,22)/RF(RFnum(22)));
    end
        
    if isempty(func)
        func = 0;
    end
end
knownMet = unique(knownMet);

end

function knownMet = findKnownMet(knownKinetics)
    knownMet = [];
    
    % dxdt1
    if all(ismember([11 12 18],knownKinetics))
        knownMet = [knownMet 1 10 14 22];
    end
    
    % dxdt2
    if all(ismember([2 3],knownKinetics))
        knownMet = [knownMet 2 8 17 21];
    end
    
    % dxdt3
    if all(ismember([15 16],knownKinetics))
        knownMet = [knownMet 3 7 9 16 22];
    end
    
    % dxdt4    
    if all(ismember([12 13],knownKinetics))
        knownMet = [knownMet 1 4 15 22];
    end
    
    % dxdt5
    if all(ismember([3 5 9 10 22 23 24],knownKinetics))
        knownMet = [knownMet 2 5 8 11 12 18 20 21];
    end
    
    % dxdt6
    if all(ismember([6 7 8],knownKinetics))
        knownMet = [knownMet 6 7 9 18 19 22];
    end
    
    % dxdt7
    if all(ismember([6 7 15],knownKinetics))
        knownMet = [knownMet 6 7 9 19 22];
    end
    
    % dxdt8
    if all(ismember([3 4 22],knownKinetics))
        knownMet = [knownMet 2 8 12 21];
    end
    
    % dxdt9
    if all(ismember([8 12 15],knownKinetics))
        knownMet = [knownMet 1 6 7 9 18 22];
    end
    
    % dxdt10
    if all(ismember([8 20],knownKinetics)) && sum(ismember([14 17 19],knownKinetics)) > 0
        knownMet = [knownMet 1 10 13];
    end
    
    % dxdt11
    if all(ismember([9 10],knownKinetics))
        knownMet = [knownMet 5 11 18 21];
    end
    
    % dxdt12
    if all(ismember([4 5],knownKinetics))
        knownMet = [knownMet 8 12 20 21];
    end
    
    % dxdt13
    if all(ismember([20 21],knownKinetics))
        knownMet = [knownMet 10 13];
    end
    
    % dxdt14
    if all(ismember([10 11],knownKinetics))
        knownMet = [knownMet 5 11 14];
    end
    
    % dxdt15
    if all(ismember([13],knownKinetics)) && sum(ismember([14 17 19],knownKinetics)) > 0
        knownMet = [knownMet 4 15];
    end
    
    % dxdt16
    if all(ismember([16],knownKinetics)) && sum(ismember([14 17 19],knownKinetics)) > 0
        knownMet = [knownMet 3 16];
    end
    
    % dxdt17
    if all(ismember([1 2],knownKinetics))
        knownMet = [knownMet 2 8 17];
    end
    
    % dxdt18
    if all(ismember([8 9],knownKinetics))
        knownMet = [knownMet 5 6 9 11 18 21 22];
    end
    
    % dxdt19
    if all(ismember([5 6],knownKinetics))
        knownMet = [knownMet 6 7 12 19 20 21];
    end
    
    % dxdt20
    if all(ismember([24],knownKinetics))
        knownMet = [knownMet 10 13 20];
    end
    
    % dxdt21
    if all(ismember([3 5 9 10 22 23 24],knownKinetics))
        knownMet = [knownMet 2 5 8 11 12 18 20 21];
    end
    
    % dxdt22
    if all(ismember([8 12 15],knownKinetics))
        knownMet = [knownMet 1 6 7 9 18 22];
    end
    knownMet = unique(knownMet);
end