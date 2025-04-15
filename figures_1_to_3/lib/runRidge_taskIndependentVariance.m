function runRidge_taskIndependentVariance(cPath,animal,rec,glmPath, fileprefix)
%for calculation of task independent and task aligned variance

NFOLDS = 10;

%get the design matrices
[regLabelsAll,regIdxAll,RA,regZeroFramesAll,zeromeanVcA,U,usedTrialsA] = ridgeModel_returnDesignMatrix_simplified(cPath,animal,rec,glmPath,'attentive',[],false,true);
[~,~,RB,~,zeromeanVcB,~,usedTrialsB] = ridgeModel_returnDesignMatrix_simplified(cPath,animal,rec,glmPath,'biased',[],false,true);

if isempty(RA) || isempty(RB) %skip sessions with too few trials
    return
end
% create some label groups
taskvarlabels = {'time', 'Choice','reward','handleSound','lfirstTacStim','lTacStim','rfirstTacStim','rTacStim','lfirstAudStim','lAudStim','rfirstAudStim','rAudStim','prevReward','prevChoice','nextChoice','water'};
allmotorlabels = {{'lGrab','lGrabRel','rGrab','rGrabRel'},{'lLick','rLick'},{'piezoAnalog','piezoDigital','piezoMoveAnalog','piezoMoveDigital','piezoMoveHiDigital'},{'whiskAnalog','whiskDigital','whiskHiDigital'},{'noseAnalog','noseDigital','noseHiDigital'},{'fastPupilAnalog','fastPupilDigital','fastPupilHiDigital','slowPupil'},{'faceAnalog','faceDigital','faceHiDigital'},{'bodyAnalog','bodyDigital','bodyHiDigital'},{'Move'},{'bhvVideo'}};
% make sure theyre in the right order and are included in the model
taskvarlabels = regLabelsAll(sort(find(ismember(regLabelsAll,taskvarlabels))));


shuffleLabels = regLabelsAll(~ismember(regLabelsAll, taskvarlabels)); %shuffle all but task vars

%run engaged trials with regressor rejection - task variable model
R = shuffleDesignMatrix(regLabelsAll,regIdxAll,RA,shuffleLabels);
[R, regLabels, regIdx, regZeroFrames, rejIdx] = ridgeModel_rejectEmptyRegressors(R,regLabelsAll,regIdxAll,regZeroFramesAll);
%% [R, regLabels, regIdx, regZeroFrames, rejIdx] = ridgeModel_rejectDeficientRegressors(R,regLabels,regIdx,regZeroFrames);
[Vm, betas, lambdas, cMap, cMovie] = ridgeModel_crossValidate(R,U,zeromeanVcA,75,NFOLDS);
ridgeModel_saveResults(cPath,animal,rec, [fileprefix 'taskA'], Vm, zeromeanVcA, U, R, betas, lambdas, cMap, cMovie, regLabels, regIdx, rejIdx, regZeroFrames, usedTrialsA);

%run disengaged trials with regressor rejection - task variable model
R = shuffleDesignMatrix(regLabelsAll,regIdxAll,RB,shuffleLabels);
[R, regLabels, regIdx, regZeroFrames, rejIdx] = ridgeModel_rejectEmptyRegressors(R,regLabelsAll,regIdxAll,regZeroFramesAll);
%%% [R, regLabels, regIdx, regZeroFrames, rejIdx] = ridgeModel_rejectDeficientRegressors(R,regLabels,regIdx,regZeroFrames);
[Vm, betas, lambdas, cMap, cMovie] = ridgeModel_crossValidate(R,U,zeromeanVcB,75,NFOLDS);
ridgeModel_saveResults(cPath,animal,rec, [fileprefix 'taskB'], Vm, zeromeanVcB, U, R, betas, lambdas, cMap, cMovie, regLabels, regIdx, rejIdx, regZeroFrames, usedTrialsB);


for i = 1:length(allmotorlabels)
    % run task models with movement variables added in
    shuffleLabels = regLabelsAll(~ismember(regLabelsAll, [taskvarlabels, allmotorlabels{i}]));

    R = shuffleDesignMatrix(regLabelsAll,regIdxAll,RA,shuffleLabels);
    [R, regLabels, regIdx, regZeroFrames, rejIdx] = ridgeModel_rejectEmptyRegressors(R,regLabelsAll,regIdxAll,regZeroFramesAll);
    %[R, regLabels, regIdx, regZeroFrames, rejIdx] = ridgeModel_rejectDeficientRegressors(R,regLabels,regIdx,regZeroFrames);
    [Vm, betas, lambdas, cMap, cMovie] = ridgeModel_crossValidate(R,U,zeromeanVcA,75,NFOLDS);
    ridgeModel_saveResults(cPath,animal,rec, [fileprefix strjoin(allmotorlabels{i}, '_') 'plustaskA'], Vm, zeromeanVcA, U, R, betas, lambdas, cMap, cMovie, regLabels, regIdx, rejIdx, regZeroFrames, usedTrialsA);

    R = shuffleDesignMatrix(regLabelsAll,regIdxAll,RB,shuffleLabels);
    [R, regLabels, regIdx, regZeroFrames, rejIdx] = ridgeModel_rejectEmptyRegressors(R,regLabelsAll,regIdxAll,regZeroFramesAll);
    %[R, regLabels, regIdx, regZeroFrames, rejIdx] = ridgeModel_rejectDeficientRegressors(R,regLabels,regIdx,regZeroFrames);
    [Vm, betas, lambdas, cMap, cMovie] = ridgeModel_crossValidate(R,U,zeromeanVcB,75,NFOLDS);
    ridgeModel_saveResults(cPath,animal,rec, [fileprefix strjoin(allmotorlabels{i}, '_') 'plustaskB'], Vm, zeromeanVcB, U, R, betas, lambdas, cMap, cMovie, regLabels, regIdx, rejIdx, regZeroFrames, usedTrialsB);

    % single variable model (for task dependent calculation)
    shuffleLabels = regLabelsAll(~ismember(regLabelsAll, allmotorlabels{i}));

    R = shuffleDesignMatrix(regLabelsAll,regIdxAll,RA,shuffleLabels);
    [R, regLabels, regIdx, regZeroFrames, rejIdx] = ridgeModel_rejectEmptyRegressors(R,regLabelsAll,regIdxAll,regZeroFramesAll);
    %[R, regLabels, regIdx, regZeroFrames, rejIdx] = ridgeModel_rejectDeficientRegressors(R,regLabels,regIdx,regZeroFrames);
    [Vm, betas, lambdas, cMap, cMovie] = ridgeModel_crossValidate(R,U,zeromeanVcA,75,NFOLDS);
    ridgeModel_saveResults(cPath,animal,rec, [fileprefix strjoin(allmotorlabels{i}, '_') 'A'], Vm, zeromeanVcA, U, R, betas, lambdas, cMap, cMovie, regLabels, regIdx, rejIdx, regZeroFrames, usedTrialsA);

    R = shuffleDesignMatrix(regLabelsAll,regIdxAll,RB,shuffleLabels);
    [R, regLabels, regIdx, regZeroFrames, rejIdx] = ridgeModel_rejectEmptyRegressors(R,regLabelsAll,regIdxAll,regZeroFramesAll);
    %[R, regLabels, regIdx, regZeroFrames, rejIdx] = ridgeModel_rejectDeficientRegressors(R,regLabels,regIdx,regZeroFrames);
    [Vm, betas, lambdas, cMap, cMovie] = ridgeModel_crossValidate(R,U,zeromeanVcB,75,NFOLDS);
    ridgeModel_saveResults(cPath,animal,rec, [fileprefix strjoin(allmotorlabels{i}, '_') 'B'], Vm, zeromeanVcB, U, R, betas, lambdas, cMap, cMovie, regLabels, regIdx, rejIdx, regZeroFrames, usedTrialsB);
end
end