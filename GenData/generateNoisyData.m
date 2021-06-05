function generateNoisyData(cleanDataName,nT,cov,rep)

% 1) Load in cleanDataName
cleanData = load(cleanDataName);

% Second '_' char should be the one at 'k-%02d_hiRes'
idxEnd = strfind(cleanDataName,'_');
saveName = sprintf('%s_nT-%03d_cov-%02d_rep-%03d.mat',cleanDataName(1:idxEnd(end)-1),nT,cov,rep);

% 2) Interpolate to lo-res nT sampling
timeVec = linspace(cleanData.tStart,cleanData.tEnd,nT+1);
fluxTimeVec = timeVec(1:end-1)+0.5*diff(timeVec(1:2));
% fluxTimeVec = timeVec;

loResConcMatrix = interp1(cleanData.timeVec,cleanData.concMatrix,timeVec,'linear','extrap');
loResFluxMatrix = interp1(cleanData.fluxTimeVec,cleanData.fluxMatrix,fluxTimeVec,'linear','extrap');


loResData.tStart = cleanData.tStart;
loResData.tEnd = cleanData.tEnd;
loResData.nT = nT;
if isfield(cleanData,'paramsVec')
    loResData.paramsVec = cleanData.paramsVec;
end
loResData.x0 = cleanData.x0;
loResData.timeVec = timeVec';
loResData.fluxTimeVec = fluxTimeVec';

% 3a) Add in noise: use noiseless data + random * cov 
loResData.concMatrix = loResConcMatrix + loResConcMatrix .* (cov/100*randn(size(loResConcMatrix)));
loResData.fluxMatrix = loResFluxMatrix + loResFluxMatrix .* (cov/101*randn(size(loResFluxMatrix)));

% Guarantee we have no *negative* concentation values
loResData.concMatrix(loResData.concMatrix<0) = 0;

% First data point is always 'correct' and noiseless
loResData.concMatrix(1,:) = loResData.x0;

% v1 is fixed value, no noise added
loResData.fluxMatrix(:,1) = loResFluxMatrix(:,1);

% 3b) Save out this noisy dataset
save(saveName,'-struct','loResData');
