function ridgeModel_saveResults(cPath,animal,rec, filename, Vm, zeromeanVc, U, R, betas, lambdas, cMap, cMovie, regLabels, regIdx, rejIdx, regZeroFrames, usedTrials)
TASK = 'SpatialDisc';

if ~strcmpi(filename(end-3:end),'.mat')
    filename = [filename '.mat'];
end
savepath = [cPath filesep animal filesep TASK filesep rec filesep filename];

save(savepath, 'Vm', 'zeromeanVc', 'U', 'R', 'betas', 'lambdas', 'cMap', 'cMovie', 'regIdx', 'regLabels', 'rejIdx', 'regZeroFrames','usedTrials');
end