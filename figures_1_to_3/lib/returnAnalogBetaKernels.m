function [regLabels, betaOut, U] = returnAnalogBetaKernels(mouse,rec,modelfile, reg)
datapath = ['X:\Widefield' filesep mouse filesep 'SpatialDisc' filesep rec filesep];
try
    load([datapath modelfile]);
catch
    fprintf('\nThe encoding model file does not exist!');
    regLabels = {};
    alignmentFrame = {};
    betaOut = {};
    U = {};
    return
end

meanbetas = mean(cat(3,betas{:}),3);

for i = 1:length(regLabels)
    betaOut{i} = meanbetas(find(regIdx == i),:);
    alignmentFrames{i} = find(regZeroFrames(find(regIdx == i))); %will be empty if there is no alignment frames
end
regind = find(strcmpi(regLabels,reg));

if isempty(regind)
    fprintf('no regressor in this model')
    fprintf('\nThe encoding model file does not exist!');
    regLabels = {};
    alignmentFrame = {};
    betaOut = {};
    U = {};
    return
end

betaOut = betaOut{regind};
end