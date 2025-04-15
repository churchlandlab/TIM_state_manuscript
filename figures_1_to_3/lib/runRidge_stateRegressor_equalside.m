function runRidge_stateRegressor_equalside(cPath,animal,rec,glmPath,fileprefix)

NFOLDS = 10;

% %get the design matrices
% [regLabelsAll,regIdxAll,RA,regZeroFramesAll,zeromeanVcA,U,usedTrialsA] = ridgeModel_returnDesignMatrix_stateReg_onestate(cPath,animal,rec,glmPath,'attentive',[],false,true);
% 
% if isempty(RA)  %skip sessions with too few trials
%     return
% end
% 
% %assert(length(intersect(usedTrialsA,usedTrialsB)) == 0,'something is wrong with trial indexing')
% 
% 
% %run full model
% [R, regLabels, regIdx, regZeroFrames, rejIdx] = ridgeModel_rejectEmptyRegressors(RA,regLabelsAll,regIdxAll,regZeroFramesAll);
% %[R, regLabels, regIdx, regZeroFrames, rejIdx] = ridgeModel_rejectDeficientRegressors(R,regLabels,regIdx,regZeroFrames);
% [Vm, betas, lambdas, cMap, cMovie] = ridgeModel_crossValidate(R,U,zeromeanVcA,75,NFOLDS);
% ridgeModel_saveResults(cPath,animal,rec, [fileprefix 'full'], Vm, zeromeanVcA, U, R, betas, lambdas, cMap, cMovie, regLabels, regIdx, rejIdx, regZeroFrames, usedTrialsA);
% 
% %run with only state - shuffle everything else
% shuffleLabels = regLabelsAll(~ismember(regLabelsAll, 'state'));
% R = shuffleDesignMatrix(regLabelsAll,regIdxAll,RA,shuffleLabels);
% 
% [R, regLabels, regIdx, regZeroFrames, rejIdx] = ridgeModel_rejectEmptyRegressors(R,regLabelsAll,regIdxAll,regZeroFramesAll);
% %[R, regLabels, regIdx, regZeroFrames, rejIdx] = ridgeModel_rejectDeficientRegressors(R,regLabels,regIdx,regZeroFrames);
% [Vm, betas, lambdas, cMap, cMovie] = ridgeModel_crossValidate(R,U,zeromeanVcA,75,NFOLDS);
% ridgeModel_saveResults(cPath,animal,rec, [fileprefix 'singlevar_state'], Vm, zeromeanVcA, U, R, betas, lambdas, cMap, cMovie, regLabels, regIdx, rejIdx, regZeroFrames, usedTrialsA);
% 
% %run without state
% shuffleLabels = regLabelsAll(ismember(regLabelsAll, 'state'));
% R = shuffleDesignMatrix(regLabelsAll,regIdxAll,RA,shuffleLabels);
% 
% [R, regLabels, regIdx, regZeroFrames, rejIdx] = ridgeModel_rejectEmptyRegressors(R,regLabelsAll,regIdxAll,regZeroFramesAll);
% %[R, regLabels, regIdx, regZeroFrames, rejIdx] = ridgeModel_rejectDeficientRegressors(R,regLabels,regIdx,regZeroFrames);
% [Vm, betas, lambdas, cMap, cMovie] = ridgeModel_crossValidate(R,U,zeromeanVcA,75,NFOLDS);
% ridgeModel_saveResults(cPath,animal,rec, [fileprefix 'no_state'], Vm, zeromeanVcA, U, R, betas, lambdas, cMap, cMovie, regLabels, regIdx, rejIdx, regZeroFrames, usedTrialsA);

%%  now run the individual states
[regLabelsAll,regIdxAll,RA,regZeroFramesAll,zeromeanVcA,U,usedTrialsA] = ridgeModel_returnDesignMatrix_simplified_onestate(cPath,animal,rec,glmPath,'attentive',[],false,true);
[regLabelsAll,regIdxAll,RB,regZeroFramesAll,zeromeanVcB,U,usedTrialsB] = ridgeModel_returnDesignMatrix_simplified_onestate(cPath,animal,rec,glmPath,'biased',[],false,true);

spontmotorlabels = {'piezoAnalog','piezoDigital','piezoMoveAnalog','piezoMoveDigital','piezoMoveHiDigital','whiskAnalog','whiskDigital','whiskHiDigital','noseAnalog','noseDigital','noseHiDigital','fastPupilAnalog','fastPupilDigital','fastPupilHiDigital','slowPupil','faceAnalog','faceDigital','faceHiDigital','bodyAnalog','bodyDigital','bodyHiDigital','Move','bhvVideo'};

if isempty(RA) || isempty(RB) %skip sessions with too few trials
    return
end

% %run engaged trials with regressor rejection
% [R, regLabels, regIdx, regZeroFrames, rejIdx] = ridgeModel_rejectEmptyRegressors(RA,regLabelsAll,regIdxAll,regZeroFramesAll);
% %[R, regLabels, regIdx, regZeroFrames, rejIdx] = ridgeModel_rejectDeficientRegressors(R,regLabels,regIdx,regZeroFrames);
% [Vm, betas, lambdas, cMap, cMovie] = ridgeModel_crossValidate(R,U,zeromeanVcA,75,NFOLDS);
% ridgeModel_saveResults(cPath,animal,rec, [fileprefix 'fullA'], Vm, zeromeanVcA, U, R, betas, lambdas, cMap, cMovie, regLabels, regIdx, rejIdx, regZeroFrames, usedTrialsA);

%run disengaged trials with regressor rejection
[R, regLabels, regIdx, regZeroFrames, rejIdx] = ridgeModel_rejectEmptyRegressors(RB,regLabelsAll,regIdxAll,regZeroFramesAll);
%[R, regLabels, regIdx, regZeroFrames, rejIdx] = ridgeModel_rejectDeficientRegressors(R,regLabels,regIdx,regZeroFrames);
[Vm, betas, lambdas, cMap, cMovie] = ridgeModel_crossValidate(R,U,zeromeanVcB,75,NFOLDS);
ridgeModel_saveResults(cPath,animal,rec, [fileprefix 'fullB'], Vm, zeromeanVcB, U, R, betas, lambdas, cMap, cMovie, regLabels, regIdx, rejIdx, regZeroFrames, usedTrialsB);

%run with only spontaneous motor variables
shuffleLabels = regLabelsAll(~ismember(regLabelsAll, spontmotorlabels));

% R = shuffleDesignMatrix(regLabelsAll,regIdxAll,RA,shuffleLabels);
% [R, regLabels, regIdx, regZeroFrames, rejIdx] = ridgeModel_rejectEmptyRegressors(R,regLabelsAll,regIdxAll,regZeroFramesAll);
% %[R, regLabels, regIdx, regZeroFrames, rejIdx] = ridgeModel_rejectDeficientRegressors(R,regLabels,regIdx,regZeroFrames);
% [Vm, betas, lambdas, cMap, cMovie] = ridgeModel_crossValidate(R,U,zeromeanVcA,75,NFOLDS);
% ridgeModel_saveResults(cPath,animal,rec, [fileprefix 'spontA'], Vm, zeromeanVcA, U, R, betas, lambdas, cMap, cMovie, regLabels, regIdx, rejIdx, regZeroFrames, usedTrialsA);

R = shuffleDesignMatrix(regLabelsAll,regIdxAll,RB,shuffleLabels);
[R, regLabels, regIdx, regZeroFrames, rejIdx] = ridgeModel_rejectEmptyRegressors(R,regLabelsAll,regIdxAll,regZeroFramesAll);
%[R, regLabels, regIdx, regZeroFrames, rejIdx] = ridgeModel_rejectDeficientRegressors(R,regLabels,regIdx,regZeroFrames);
[Vm, betas, lambdas, cMap, cMovie] = ridgeModel_crossValidate(R,U,zeromeanVcB,75,NFOLDS);
ridgeModel_saveResults(cPath,animal,rec, [fileprefix 'spontB'], Vm, zeromeanVcB, U, R, betas, lambdas, cMap, cMovie, regLabels, regIdx, rejIdx, regZeroFrames, usedTrialsB);


end