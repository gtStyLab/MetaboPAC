function [absolute_concMatrix predicted_responseFactors] = MetaboPAC_ecoli_log(rep,percKnownKinetics)

% Set up initial information about system
modelInfo_chass
load('chassV_k-01_hiRes.mat');

modelInfo.fixedFluxes = (modelInfo.vBounds(:,1)==modelInfo.vBounds(:,2));

numMetabs = size(modelInfo.S,1);
numFlux = size(modelInfo.S,2);

numExtraMets = numMetabs - size(concMatrix,2);
if numExtraMets > 0
    concMatrix(:,end+1:end+numExtraMets) = zeros(size(concMatrix,1),numExtraMets);
end

% Find Mass Action reactions
[row_MA col_MA] = find(modelInfo.S < 0);
MA_reactions = [row_MA col_MA];

% Load relative abundances
load('chass_trueRF_log');
trueRF = trueRF(rep,:);

% Calculate relative abundances
relative_concMatrix = concMatrix.*trueRF;

% Perform kinetic equations approach
numKnownKinetics = round(percKnownKinetics/100*numFlux);
knownKinetics = sort(randperm(numFlux,numKnownKinetics));
for i = 1:48
    [KEapproach_results(i,:) knownMet] = solveRF_KEapproach_chass(relative_concMatrix,timeVec,knownKinetics);
end
if knownMet ~= 0
    RF_kinetics = median(KEapproach_results);
else
    RF_kinetics = [];
    knownMet = [];
end

% Set up information for optimization approach
maxRandVal = 1000;
minRandVal = 1;
lb = minRandVal*ones(1,numMetabs-length(knownMet));
ub = maxRandVal*ones(1,numMetabs-length(knownMet));
trueRF_lower_mag = 0;
trueRF_upper_mag = 3;
options.ConstraintTolerance	= 1e10;
options.MaxFunctionEvaluations = 5000;

% Calculate relative pooling fluxes by calculating change in relative
% abundances and dividing by change in time
for i = 2:size(relative_concMatrix,1)
    relative_Vpool(i-1,:) = (relative_concMatrix(i,:)-relative_concMatrix(i-1,:))./(timeVec(i)-timeVec(i-1));
end

% Perform optimization approach
fval = 1500*ones(48,1);
RF_opt = ones(48,numMetabs);
if length(knownMet) ~= numMetabs
    for i = 1:48
        % Set initial seed for optimizier
        x0 = 10.^(trueRF_lower_mag + (trueRF_upper_mag-trueRF_lower_mag)*rand(1,numMetabs-length(knownMet)));
        ga_options = optimoptions('ga','MaxTime',30,'InitialPopulationMatrix',x0);

        % Perform genetic algorithm optimization
        [optimalRF fval(i,1)] = ga(@(testRF) calcPenalty(testRF,modelInfo,relative_Vpool,numMetabs,numFlux,MA_reactions,relative_concMatrix,timeVec,fluxTimeVec,RF_kinetics,knownMet,knownKinetics),numMetabs-length(knownMet),[],[],[],[],lb,ub,[],ga_options);
        
        RF_temp = zeros(1,numMetabs);
        RF_temp(1,knownMet) = RF_kinetics;
        RF_temp(1,setdiff(1:numMetabs,knownMet)) = optimalRF;

        RF_opt(i,:) = RF_temp;
    end
else
    % If optimization approach not performed, use results from kinetic
    % equations approach
    RF_opt = KEapproach_results;
end

predicted_responseFactors = median(RF_opt);
absolute_concMatrix = relative_concMatrix./predicted_responseFactors;

save(sprintf('ecoli_MetaboPAC_log_percKnownKinetics-%03d_rep-%03d.mat',percKnownKinetics,rep));
end



function penalty = calcPenalty(testRF,modelInfo,relative_Vpool,numMetabs,numFlux,MA_reactions,relative_concMatrix,timeVec,fluxTimeVec,RF_kinetics,knownMet,knownKinetics)
    % Consolidate response factor values
    RF = zeros(1,numMetabs);
    RF(1,knownMet) = RF_kinetics;
    RF(1,setdiff(1:numMetabs,knownMet)) = testRF;

    % Calculate change in inferred absolute concentrations over time (i.e.
    % pooling fluxes, Vpool)
    test_Vpool(:,1:numMetabs) = relative_Vpool(:,1:numMetabs)./RF;
    
    % Calculate inferred absolute concentrations
    absolute_concMatrix(:,1:numMetabs) = relative_concMatrix(:,1:numMetabs)./RF;
    
    % Set up flux matrix before inferring fluxes
    fluxMatrixTemp = nan(size(test_Vpool,1),numFlux+numMetabs);
    fluxMatrixTemp(:,numFlux+1:end) = test_Vpool;
    fluxMatrixTemp(:,modelInfo.fixedFluxes) = ones(size(fluxMatrixTemp(:,modelInfo.fixedFluxes),1),1)*modelInfo.vBounds(modelInfo.fixedFluxes,1)';
    newFluxMatrixTemp = calcFluxWithKinetics_chass(absolute_concMatrix,timeVec,knownKinetics);
    modelInfo.fixedFluxes = ~isnan(newFluxMatrixTemp(1,1:numFlux))';
    
    % Infer fluxes using pinv
    Vcalc = calcFluxesViaPinv(newFluxMatrixTemp,modelInfo.S,modelInfo.fixedFluxes);    
    
    % Calculate penalties
    if sum(sum(Vcalc))~=0 && ~isnan(sum(sum(Vcalc)))
        % Calculate pooling flux penalty
        for i = 1:size(Vcalc,1)
            Vcalc_pool(i,:) = modelInfo.S*Vcalc(i,1:numFlux)';
        end

        % Calculate mass balance penalty
        massbalance_penalty = sqrt(sum(sum((Vcalc_pool(2:end,:) - Vcalc(2:end,numFlux+1:end)).^2)));
        %rmse = sqrt(mean((Vcalc_pool(2:end,:) - Vcalc(2:end,numFlux+1:end)).^2,'all'));
        %nrmse = sqrt(mean(((Vcalc_pool(2:end,:) - Vcalc(2:end,numFlux+1:end))./range(Vcalc(2:end,numFlux+1:end))).^2,'all'));
        
        % Calculate max. concentration penalty
        if any(any(absolute_concMatrix(:,1:numMetabs-2) > 50))
            conc_penalty = max(max(absolute_concMatrix(:,1:numMetabs-2)));
        else
            conc_penalty = 0;
        end
        

        % Calculate correlation penalty for reactions controlled by one
        % metabolite
        oneContMetCorr_penalty = pen_oneContMetCorr_chass(absolute_concMatrix,Vcalc,timeVec,fluxTimeVec);
        
        % Calculate correlation penalty for reactions controlled by two
        % metabolites
        twoContMetCorr_penalty = pen_twoContMetCorr_chass(absolute_concMatrix,Vcalc,timeVec,fluxTimeVec);

        % Calculate curve fit penalty for reactions controlled by one
        % metabolite
        oneContMetCurveFit_penalty = pen_oneContMetCurveFit_chass(absolute_concMatrix,Vcalc,timeVec,fluxTimeVec);

        % Calculate BST fit penalty
        nT = size(absolute_concMatrix,1)-1;
        BST_penalty = pen_BSTfit_chass(nT,absolute_concMatrix,Vcalc);
        
        % Calculate total penalty
        penalty = abs(massbalance_penalty)/1000 + abs(conc_penalty) + abs(oneContMetCorr_penalty) + abs(twoContMetCorr_penalty) + abs(oneContMetCurveFit_penalty) + abs(BST_penalty);
        
        if isnan(penalty)
            penalty = 1e7;
        end
    else
        penalty = sum(1./RF)*1e7;
    end
end
