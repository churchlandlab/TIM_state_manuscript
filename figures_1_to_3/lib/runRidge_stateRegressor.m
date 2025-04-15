function runRidge_stateRegressor(cPath,animal,rec,glmPath,fileprefix)

NFOLDS = 10;

%get the design matrices
[regLabelsAll,regIdxAll,RA,regZeroFramesAll,zeromeanVcA,U,usedTrialsA] = ridgeModel_returnDesignMatrix_stateReg(cPath,animal,rec,glmPath,'attentive',[],false,true);
if isempty(RA) %skip sessions with too few trials
    return
end

%assert(length(intersect(usedTrialsA,usedTrialsB)) == 0,'something is wrong with trial indexing')


%run full model
[R, regLabels, regIdx, regZeroFrames, rejIdx] = ridgeModel_rejectEmptyRegressors(RA,regLabelsAll,regIdxAll,regZeroFramesAll);
%[R, regLabels, regIdx, regZeroFrames, rejIdx] = ridgeModel_rejectDeficientRegressors(R,regLabels,regIdx,regZeroFrames);
[Vm, betas, lambdas, cMap, cMovie] = ridgeModel_crossValidate(R,U,zeromeanVcA,75,NFOLDS);
ridgeModel_saveResults(cPath,animal,rec, [fileprefix 'full'], Vm, zeromeanVcA, U, R, betas, lambdas, cMap, cMovie, regLabels, regIdx, rejIdx, regZeroFrames, usedTrialsA);

%run with only state - shuffle everything else
shuffleLabels = regLabelsAll(~ismember(regLabelsAll, 'state'));
R = shuffleDesignMatrix(regLabelsAll,regIdxAll,RA,shuffleLabels);

[R, regLabels, regIdx, regZeroFrames, rejIdx] = ridgeModel_rejectEmptyRegressors(R,regLabelsAll,regIdxAll,regZeroFramesAll);
%[R, regLabels, regIdx, regZeroFrames, rejIdx] = ridgeModel_rejectDeficientRegressors(R,regLabels,regIdx,regZeroFrames);
[Vm, betas, lambdas, cMap, cMovie] = ridgeModel_crossValidate(R,U,zeromeanVcA,75,NFOLDS);
ridgeModel_saveResults(cPath,animal,rec, [fileprefix 'singlevar_state'], Vm, zeromeanVcA, U, R, betas, lambdas, cMap, cMovie, regLabels, regIdx, rejIdx, regZeroFrames, usedTrialsA);

%run without state
shuffleLabels = regLabelsAll(ismember(regLabelsAll, 'state'));
R = shuffleDesignMatrix(regLabelsAll,regIdxAll,RA,shuffleLabels);

[R, regLabels, regIdx, regZeroFrames, rejIdx] = ridgeModel_rejectEmptyRegressors(R,regLabelsAll,regIdxAll,regZeroFramesAll);
%[R, regLabels, regIdx, regZeroFrames, rejIdx] = ridgeModel_rejectDeficientRegressors(R,regLabels,regIdx,regZeroFrames);
[Vm, betas, lambdas, cMap, cMovie] = ridgeModel_crossValidate(R,U,zeromeanVcA,75,NFOLDS);
ridgeModel_saveResults(cPath,animal,rec, [fileprefix 'no_state'], Vm, zeromeanVcA, U, R, betas, lambdas, cMap, cMovie, regLabels, regIdx, rejIdx, regZeroFrames, usedTrialsA);

end