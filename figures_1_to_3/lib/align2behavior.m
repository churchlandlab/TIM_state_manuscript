function [alignVc,bhv,goodtrials] = align2behavior(cPath,Animal,Rec,trialInds,hemoCorrected)
%where trialInds is the desired trials to select for alignment
%pass trialinds based on the sessiondata file, this function will work out
%the imaging data itself and discard requested trials that dont have imaging data.
addpath('C:\Data\churchland\ridgeModel\widefield');
addpath('C:\Data\churchland\ridgeModel\smallstuff');
Paradigm = 'SpatialDisc';
cPath = [cPath filesep Animal filesep Paradigm filesep Rec filesep]; %Widefield data path
dims = 200; %number of dims of widefield SVD
%% Load data
bhvFile = dir([cPath filesep Animal '_' Paradigm '*.mat']);
load([cPath bhvFile(1).name],'SessionData'); %load behavior data

if hemoCorrected
    load([cPath 'Vc.mat'],'Vc','U','trials','bTrials'); %just need trials variable here. Vs is [dims of temporal components,frames,trials]
else
    load([cPath 'Vc.mat'],'trials','bTrials');
    %load([cPath 'hemoV.mat'], 'hemoV','U');
    %Vc = hemoV(1:200,:,:);
    load([cPath 'blueV.mat'], 'blueV','U');
    Vc = blueV(1:200,:,:);
    U = U(:,:,1:200);
    clear hemoV;
end


%%%%%%%

try
    load([cPath 'opts3.mat']);
    opts = opts3;clear opts3;
catch
    try
        load([cPath 'opts2.mat']);
    catch
        load([cPath 'opts.mat']);
    end
end
sRate = opts.frameRate;
nRequested = length(trialInds);

%% get proper trials from Vc and SessionData
% ensure there are not too many trials in Vc
ind = trials > SessionData.nTrials;
trials(ind) = [];
bTrials(ind) = []; %btrials gives the indices of bpod trial numbers that have imaging data
Vc(:,:,ind) = [];


temp = trialInds(ismember(trialInds,bTrials)); %use only trial indices to grab that are in the Vc dataset. contains only bpod trial numbers that have imaging data.
bhv = selectBehaviorTrials(SessionData,temp); %this grabs trials that are in the requested indices AND the Vc dataset
bTrials = find(ismember(bTrials,trialInds)); %gets the indices of trialInds that have their bpod trial number contained in bTrials
Vc = Vc(:,:,bTrials);
nTrials = length(bTrials);
goodtrials = temp;

fprintf('Requested %i trials. Using %i trials (Vc dataset usually has a few missing trials).\n',nRequested,nTrials);

[~,~,a] = size(Vc);
b = length(bhv.CorrectSide);
if a~=b
    error('Size of data is not consistent!!!')
end
%segIdx = [2 0.75 1.25 0.75 1]; %what simon sent me.
if sum(ismember(Animal,'mSM')) == 3 %mSM Mice
    segIdx = [1 0.5 1.00 0.75 .75]; %[baseline, handle, stim, delay, response] maximal duration of each segment in seconds, use this for EMX mice
    segIdx = [1 0.5 1.00 0.4 .75] %testing, used for decoder to keep shuffled decoder distribution to chance
elseif sum(ismember(Animal,'Fez')) == 3 %Fez mice
    segIdx = [1 0.5 1.00 0.4 .75] %testing, used for decoder to keep shuffled decoder distribution to chance
elseif sum(ismember(Animal,'CSP')) == 3 %CSP Mice
    segIdx = [1 0.2 .5 0.15 .75]; %[baseline, handle, stim, delay, response] maximal duration of each segment in seconds, use this for CSP mice
    segIdx = [1 0.5 1.00 0.4 .75] %testing
end

segFrames = cumsum(floor(segIdx * sRate));
alignVc.all = rateDisc_getBhvRealignment(Vc, bhv, segFrames, opts);
alignVc.segIdx = segIdx;
alignVc.segFrames = segFrames;
alignVc.transParams = opts.transParams;
alignVc.U = U;
alignVc.fs = opts.frameRate;
end
