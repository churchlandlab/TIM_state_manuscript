function [regLabels,regIdx,fullR,regZeroFrames,zeromeanVc,U,usedTrials] = ridgeModel_returnDesignMatrix_stateReg(cPath,Animal,Rec,glmPath,desiredstate,dType,randomizeStates,orthogvideotoinstructed)
%By Max Melin. Trains the ridge regression model described in Musall 2019
%on that same dataset. In this code, no state regressors are added. Rather,
%trials are split up by state (predicted by the Ashwood GLM-HMM) and used to
%train seperate encoding models for each behavioral state. A model is built
%for each of the following states: ___,___,___.

% desiredstate is 1,2 or 3 (attentive, lbias, or r bias) can pass an array
% of multiple

if ~strcmpi(cPath(end),filesep)
    cPath = [cPath filesep];
end

if ~exist('dType', 'var') || isempty(dType)
    dType = 'Widefield'; %default is widefield data
end

if strcmpi(dType,'twoP')
    piezoLine = 5;     % channel in the analog data that contains data from piezo sensor
    stimLine = 4;      % channel in the analog data that contains stimulus trigger.

elseif strcmpi(dType,'Widefield')
    piezoLine = 2;     % channel in the analog data that contains data from piezo sensor
    stimLine = 6;      % channel in the analog data that contains stimulus trigger.
end

Paradigm = 'SpatialDisc';
cPath = [cPath Animal filesep Paradigm filesep Rec filesep]; %Widefield data path
sPath = ['\\grid-hs\churchland_nlsas_data\data\BpodImager\Animals\' Animal filesep Paradigm filesep Rec filesep]; %server data path. not used on hpc.
% sPath = ['/sonas-hs/churchland/hpc/home/space_managed_data/BpodImager/Animals/' Animal filesep Paradigm filesep Rec filesep]; %server data path. not used on hpc.

%% general variables
try
    load([cPath 'opts.mat'], 'opts');         % get some options from imaging data
catch
    fprintf('no opts file\n')
    regIdx=[];
    regLabels = [];
    fullR = [];
    regZeroFrames = [];
    zeromeanVc = [];
    U = [];
    usedTrials = [];
    return
end

sRate = opts.frameRate;           % Sampling rate of imaging in Hz
if sRate == 30
    vcFile = 'rsVc.mat'; %use downsampled data here and change sampling frequency
    sRate = 15;
else
    vcFile = 'Vc.mat';
end

preStimDur = floor(2 * sRate) / sRate; % Duration of trial before lever grab in seconds
postStimDur = floor(3 * sRate) / sRate; % Duration of trial after lever grab onset in seconds
frames = round((preStimDur + postStimDur) * sRate); %nr of frames per trial
trialDur = (frames * (1/sRate)); %duration of trial in seconds

%other variables
mPreTime = ceil(0.5 * sRate);  % precede motor events to capture preparatory activity in frames
mPostTime = ceil(2 * sRate);   % follow motor events for mPostStim in frames
fsPostTime = ceil(3 * sRate);   % follow first stim event for fsPostStim in frames
sPostTime = ceil(2 * sRate);   % follow stim events for sPostStim in frames
motorIdx = [-(mPreTime: -1 : 1) 0 (1:mPostTime)]; %index for design matrix to cover pre- and post motor action
tapDur = 0.1;      % minimum time of lever contact, required to count as a proper grab.
leverMoveDur = 0.3; %duration of lever movement. this is used to orthogonalize video against lever movement.
leverMoveDur = ceil(leverMoveDur * sRate); %convert to frames
shVal = sRate * opts.preStim  + 1; %expected position of stimulus onset in the imaging data (s).
maxStimShift = 1 * sRate; % maximal stimulus onset after handle grab. (default is 1s - this means that stimulus onset should not be more than 1s after handle grab. Cognitive regressors will have up to 1s of baseline because of stim jitter.)
bhvDimCnt = 200;    % number of dimensions from behavioral videos that are used as regressors.
gaussShift = 1;     % inter-frame interval between regressors. Will use only every 'gaussShift' regressor and convolve with gaussian of according FHWM to reduce total number of used regressors.
dims = 200; %number of dims in Vc that are used in the model


%% load data
bhvFile = dir([cPath filesep Animal '_' Paradigm '*.mat']);
load([cPath bhvFile(1).name],'SessionData'); %load behavior data
load(glmPath,'mouse','masks','glmhmm_params','posterior_probs','model_training_sessions','state_label_indices'); %load latent states and model info
SessionData.TrialStartTime = SessionData.TrialStartTime * 86400; %convert trailstart timestamps to seconds
nochoice = isnan(SessionData.ResponseSide); %trials without choice. used for interpolation of latent state on NaN choice trials (GLMHMM doesn't predict for these trials)

model_training_sessions = cellstr(model_training_sessions); %convert to cell
mouse = cellstr(mouse);
sessionidx = find(strcmp(Rec,model_training_sessions) & strcmp(Animal,mouse)); %find the index of the session we want to pull latent states for
postprob_nonan = posterior_probs{sessionidx}; %get latent states for desired session
mask = masks{sessionidx};
counterind = 1;
%this should probably be moved to the python side
postprob_withnan = NaN(1,size(postprob_nonan,2));
for i = 1:length(mask) %this for loop adds nan's to the latent state array. The nans will ultimatel get discarded later since the encoding model doesn't use trials without choice.
    if mask(i) == 1 %if not masked
        postprob_withnan(i,:) = postprob_nonan(counterind,:); %just put the probabilities into the new array
        counterind = counterind + 1;
    else %if no choice was made
        postprob_withnan(i,:) = NaN; %insert a NaN to new array
    end
end
postprobs = postprob_withnan;
postprobs = postprobs(:,state_label_indices); %permute the states to the correct indices
clear postprob_withnan postprob_nonan;
if randomizeStates
    postprobs = postprobs(randperm(size(postprobs,1)),:);
end
glmweights = squeeze(glmhmm_params.observations.Wk);
glmweights = glmweights(state_label_indices,:);
fprintf('\nGLM weights:\n');
disp(glmweights);

% figure;
% plot(glmweights','LineWidth',2);%plot states to make sure predicted labels are reasonable
% y = 5;
% line([1,2],[0,0],'LineStyle','--','Color','black');
% title('Latent state weights - Verify that these labels look correct');
% legend('Attentive','Lbias','Rbias');
% drawnow;
%
% figure;
% plot(postprobs,'LineWidth',2);%plot states to make sure predicted labels are reasonable
% y = 5;
% title('Posterior probabilities');
% legend('Attentive','Lbias','Rbias');
% drawnow;


if exist([cPath vcFile],'file') ~= 2 %check if data file exists and get from server otherwise
    copyfile([sPath vcFile],[cPath vcFile]);
    copyfile([sPath 'mask.mat'],[cPath 'mask.mat']);
    bhvFile = dir([sPath filesep Animal '_' Paradigm '*.mat']);
    copyfile([sPath bhvFile.name],[cPath bhvFile.name]);
end
load([cPath vcFile], 'Vc', 'U', 'trials', 'bTrials')
Vc = Vc(1:dims,:,:);
U = U(:,:,1:dims);

% ensure there are not too many trials in Vc
ind = trials > SessionData.nTrials;
trials(ind) = [];
bTrials(ind) = [];
Vc(:,:,ind) = [];

bhv = selectBehaviorTrials(SessionData,bTrials); %only use completed trials that are in the Vc dataset
postprobs = postprobs(bTrials,:); %trim to sessions with good imaging


%%%%%%%%%%%%%%%%% this just gets the state with the highest probability
[~,inds] = max(postprobs,[],2);
attentiveind = (inds == 1)'; %get attentive trials
biasind = (inds ~= 1)';
useIdx = ~isnan(bhv.ResponseSide); %only use performed trials
%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%% this gets states greater than .8
% binary = postprobs > .8; %get indices with P(state) > .8
% for i = 1:size(binary,1)
%     temp = binary(i,:) == 1;
%     if sum(temp) == 0
%         state1hot(i) = 0; % if no states are P > .8, set to 0 to indicate no state
%     else
%         state1hot(i) = find(temp);
%     end
% end
% attentiveind = state1hot == 1';
% biasind = state1hot == 2 | state1hot == 3';
% useIdx = ~isnan(bhv.ResponseSide); %only use performed trials
% useIdx = useIdx & (attentiveind | biasind); %and only use trials where P of ANY state was > .8
% % %%%%%%%%%%%%%%%%%%%

%equalize state and reward
choiceIdx = rateDisc_equalizeTrials(useIdx, attentiveind, bhv.Rewarded, inf, true);



trials = trials(choiceIdx);
bTrials = bTrials(choiceIdx);
Vc = Vc(:,:,choiceIdx);
bhv = selectBehaviorTrials(SessionData,bTrials); %only use completed trials that are in the Vc dataset
trialCnt = length(bTrials);
usedTrials = bTrials;
states = attentiveind(choiceIdx);
%mintrials = 30;
mintrials = 50;

if trialCnt < mintrials
    fprintf('\there are less than %i trials, aborting for this session.\n',mintrials);
    regIdx=[];
    regLabels = [];
    fullR = [];
    regZeroFrames = [];
    zeromeanVc = [];
    U = [];
    return
end
%% load behavior data
if exist([cPath 'BehaviorVideo' filesep 'SVD_CombinedSegments.mat'],'file') ~= 2 || ... %check if svd behavior exists on hdd and pull from server otherwise
        exist([cPath 'BehaviorVideo' filesep 'motionSVD_CombinedSegments.mat'],'file') ~= 2

    if ~exist([cPath 'BehaviorVideo' filesep], 'dir')
        mkdir([cPath 'BehaviorVideo' filesep]);
    end
    copyfile([sPath 'BehaviorVideo' filesep 'SVD_CombinedSegments.mat'],[cPath 'BehaviorVideo' filesep 'SVD_CombinedSegments']);
    copyfile([sPath 'BehaviorVideo' filesep 'motionSVD_CombinedSegments.mat'],[cPath 'BehaviorVideo' filesep 'motionSVD_CombinedSegments.mat']);
    copyfile([sPath 'BehaviorVideo' filesep 'FilteredPupil.mat'],[cPath 'BehaviorVideo' filesep 'FilteredPupil.mat']);
    copyfile([sPath 'BehaviorVideo' filesep 'segInd1.mat'],[cPath 'BehaviorVideo' filesep 'segInd1.mat']);
    copyfile([sPath 'BehaviorVideo' filesep 'segInd2.mat'],[cPath 'BehaviorVideo' filesep 'segInd2.mat']);

    movFiles = dir([sPath 'BehaviorVideo' filesep '*Video_*1.mj2']);
    copyfile([sPath 'BehaviorVideo' filesep movFiles(1).name],[cPath 'BehaviorVideo' filesep movFiles(1).name]);
    movFiles = dir([sPath 'BehaviorVideo' filesep '*Video_*2.mj2']);
    copyfile([sPath 'BehaviorVideo' filesep movFiles(1).name],[cPath 'BehaviorVideo' filesep movFiles(1).name]);

    svdFiles = dir([sPath 'BehaviorVideo' filesep '*SVD*-Seg*']);
    for iFiles = 1:length(svdFiles)
        copyfile([sPath 'BehaviorVideo' filesep svdFiles(iFiles).name],[cPath 'BehaviorVideo' filesep svdFiles(iFiles).name]);
    end
end

load([cPath 'BehaviorVideo' filesep 'SVD_CombinedSegments.mat'],'vidV'); %load behavior video data
V1 = vidV(:,1:bhvDimCnt); %behavioral video regressors
load([cPath 'BehaviorVideo' filesep 'motionSVD_CombinedSegments.mat'],'vidV'); %load abs motion video data
V2 = vidV(:,1:bhvDimCnt); % motion regressors

% check options that were used for dimensionality reduction and ensure that imaging and video data trials are equal length
load([cPath 'BehaviorVideo' filesep 'bhvOpts.mat'], 'bhvOpts'); %load abs motion video data
bhvRate = bhvOpts.targRate; %framerate of face camera
if (bhvOpts.preStimDur + bhvOpts.postStimDur) > (opts.preStim + opts.postStim) %if behavioral video trials are longer than imaging data
    V1 = reshape(V1, [], SessionData.nTrials, bhvDimCnt);
    V2 = reshape(V2, [], SessionData.nTrials, bhvDimCnt);
    if bhvOpts.preStimDur > opts.preStim
        frameDiff = ceil((bhvOpts.preStimDur - opts.preStim) * bhvRate); %overhead in behavioral frames that can be removed.
        V1 = V1(frameDiff+1:end, :, :); %cut data to size
        V2 = V2(frameDiff+1:end, :, :);
    end
    if bhvOpts.postStimDur > opts.postStim
        frameDiff = ceil((bhvOpts.postStimDur - opts.postStim) * bhvRate); %overhead in behavioral frames that can be removed.
        V1 = V1(1 : end - frameDiff, :, :); %cut data to size
        V2 = V2(1 : end - frameDiff, :, :);
    end
    V1 = reshape(V1, [], bhvDimCnt);
    V2 = reshape(V2, [], bhvDimCnt);
end

load([cPath 'BehaviorVideo' filesep 'FilteredPupil.mat'], 'pTime', 'fPupil', 'sPupil', 'whisker', 'faceM', 'bodyM', 'nose', 'bTime'); %load pupil data
%check if timestamps from pupil data are shifted against bhv data
timeCheck1 = (SessionData.TrialStartTime(1)) - (pTime{1}(1)); %time difference between first acquired frame and onset of first trial
timeCheck2 = (SessionData.TrialStartTime(1)) - (bTime{1}(1)); %time difference between first acquired frame and onset of first trial
if (timeCheck1 > 3590 && timeCheck1 < 3610) && (timeCheck2 > 3590 && timeCheck2 < 3610) %timeshift by one hour (+- 10seconds)
    warning('Behavioral and video timestamps are shifted by 1h. Will adjust timestamps in video data accordingly.')
    for iTrials = 1 : length(pTime)
        pTime{iTrials} = pTime{iTrials} + 3600; %add one hour
        bTime{iTrials} = bTime{iTrials} + 3600; %add one hour
    end
elseif timeCheck1 > 30 || timeCheck1 < -30 || timeCheck2 > 30 || timeCheck2 < -30
    error('Something wrong with timestamps in behavior and video data. Time difference is larger as 30 seconds.')
end

if any(bTrials > length(pTime))
    warning(['There are insufficient trials in the pupil data. Rejected the last ' num2str(sum(bTrials > length(pTime))) ' trial(s)']);
    bTrials(bTrials > length(pTime)) = [];
    trialCnt = length(bTrials);
end


%% find events in BPod time - All timestamps are relative to stimulus onset event to synchronize to imaging data later
% pre-allocate vectors
lickL = cell(1,trialCnt);
lickR = cell(1,trialCnt);
leverIn = NaN(1,trialCnt);
levGrabL = cell(1,trialCnt);
levGrabR = cell(1,trialCnt);
levReleaseL = cell(1,trialCnt);
levReleaseR = cell(1,trialCnt);
water = NaN(1,trialCnt);
handleSounds = cell(1,trialCnt);

tacStimL = cell(1,trialCnt);
tacStimR = cell(1,trialCnt);
audStimL = cell(1,trialCnt);
audStimR = cell(1,trialCnt);

for iTrials = 1:trialCnt

    leverTimes = [reshape(bhv.RawEvents.Trial{iTrials}.States.WaitForAnimal1',1,[]) ...
        reshape(bhv.RawEvents.Trial{iTrials}.States.WaitForAnimal2',1,[]) ...
        reshape(bhv.RawEvents.Trial{iTrials}.States.WaitForAnimal3',1,[])];

    try
        stimGrab(iTrials) = leverTimes(find(leverTimes == bhv.RawEvents.Trial{iTrials}.States.WaitForCam(1))-1); %find start of lever state that triggered stimulus onset
        handleSounds{iTrials} = leverTimes(1:2:end) - stimGrab(iTrials); %track indicator sound when animal is grabing both handles
        stimTime(iTrials) = bhv.RawEvents.Trial{iTrials}.Events.Wire3High - stimGrab(iTrials); %time of stimulus onset - measured from soundcard
        stimEndTime(iTrials) = bhv.RawEvents.Trial{iTrials}.States.DecisionWait(1) - stimGrab(iTrials); %end of stimulus period, relative to handle grab, from simons email
    catch
        stimTime(iTrials) = NaN;
        stimGrab(iTrials) = 0;
    end

    %check for spout motion
    if isfield(bhv.RawEvents.Trial{iTrials}.States,'MoveSpout')
        spoutTime(iTrials) = bhv.RawEvents.Trial{iTrials}.States.MoveSpout(1) - stimGrab(iTrials);

        %also get time when the other spout was moved out at
        if bhv.Rewarded(iTrials)
            spoutOutTime(iTrials) = bhv.RawEvents.Trial{iTrials}.States.Reward(1) - stimGrab(iTrials);
        else
            spoutOutTime(iTrials) = bhv.RawEvents.Trial{iTrials}.States.HardPunish(1) - stimGrab(iTrials);
        end
    else
        spoutTime(iTrials) = NaN;
        spoutOutTime(iTrials) = NaN;
    end

    if isfield(bhv.RawEvents.Trial{iTrials}.Events,'Port1In') %check for licks
        lickL{iTrials} = bhv.RawEvents.Trial{iTrials}.Events.Port1In;
        lickL{iTrials}(lickL{iTrials} < bhv.RawEvents.Trial{iTrials}.States.MoveSpout(1)) = []; %dont use false licks that occured before spouts were moved in
        lickL{iTrials} = lickL{iTrials} - stimGrab(iTrials);
    end
    if isfield(bhv.RawEvents.Trial{iTrials}.Events,'Port3In') %check for right licks
        lickR{iTrials} = bhv.RawEvents.Trial{iTrials}.Events.Port3In;
        lickR{iTrials}(lickR{iTrials} < bhv.RawEvents.Trial{iTrials}.States.MoveSpout(1)) = []; %dont use false licks that occured before spouts were moved in
        lickR{iTrials} = lickR{iTrials} - stimGrab(iTrials);
    end

    % get stimulus events times
    audStimL{iTrials} = bhv.stimEvents{iTrials}{1} + stimTime(iTrials);
    audStimR{iTrials} = bhv.stimEvents{iTrials}{2} + stimTime(iTrials);
    tacStimL{iTrials} = bhv.stimEvents{iTrials}{5} + stimTime(iTrials);
    tacStimR{iTrials} = bhv.stimEvents{iTrials}{6} + stimTime(iTrials);

    leverIn(iTrials) = min(bhv.RawEvents.Trial{iTrials}.States.Reset(:)) - stimGrab(iTrials); %first reset state causes lever to move in

    if isfield(bhv.RawEvents.Trial{iTrials}.Events,'Wire2High') %check for left grabs
        levGrabL{iTrials} = bhv.RawEvents.Trial{iTrials}.Events.Wire2High - stimGrab(iTrials);
    end
    if isfield(bhv.RawEvents.Trial{iTrials}.Events,'Wire1High') %check for right grabs
        levGrabR{iTrials} = bhv.RawEvents.Trial{iTrials}.Events.Wire1High - stimGrab(iTrials);
    end

    if isfield(bhv.RawEvents.Trial{iTrials}.Events,'Wire2Low') %check for left release
        levReleaseL{iTrials} = bhv.RawEvents.Trial{iTrials}.Events.Wire2Low - stimGrab(iTrials);
    end
    if isfield(bhv.RawEvents.Trial{iTrials}.Events,'Wire1Low') %check for right release
        levReleaseR{iTrials} = bhv.RawEvents.Trial{iTrials}.Events.Wire1Low - stimGrab(iTrials);
    end

    if ~isnan(bhv.RawEvents.Trial{iTrials}.States.Reward(1)) %check for reward state
        water(iTrials) = bhv.RawEvents.Trial{iTrials}.States.Reward(1) - stimGrab(iTrials);
    end
end
maxSpoutRegs = length(min(round((preStimDur + spoutTime) * sRate)) : frames); %maximal number of required spout regressors

%% build regressors - create design matrix based on event times
%basic time regressors
timeR = logical(diag(ones(1,frames)));

lGrabR = cell(1,trialCnt);
lGrabRelR = cell(1,trialCnt);
rGrabR = cell(1,trialCnt);
rGrabRelR = cell(1,trialCnt);
lLickR = cell(1,trialCnt);
rLickR = cell(1,trialCnt);
leverInR = cell(1,trialCnt);

lfirstTacStimR = cell(1,trialCnt);
rfirstTacStimR = cell(1,trialCnt);
lfirstAudStimR = cell(1,trialCnt);
rfirstAudStimR = cell(1,trialCnt);

lTacStimR = cell(1,trialCnt);
rTacStimR = cell(1,trialCnt);
lAudStimR = cell(1,trialCnt);
rAudStimR = cell(1,trialCnt);

spoutR = cell(1,trialCnt);
spoutOutR = cell(1,trialCnt);

rewardR = cell(1,trialCnt);
ChoiceR = cell(1,trialCnt);
stateR = cell(1,trialCnt);

prevRewardR = cell(1,trialCnt);
prevChoiceR = cell(1,trialCnt);
prevStimR = cell(1,trialCnt);
nextChoiceR = cell(1,trialCnt);
repeatChoiceR = cell(1,trialCnt);

waterR = cell(1,trialCnt);
fastPupilR = cell(1,trialCnt);
slowPupilR = cell(1,trialCnt);

whiskR = cell(1,trialCnt);
noseR = cell(1,trialCnt);
piezoR = cell(1,trialCnt);
piezoMoveR = cell(1,trialCnt);
faceR = cell(1,trialCnt);
bodyR = cell(1,trialCnt);

handleSoundR = cell(1,trialCnt);

%%
tic
for iTrials = 1:trialCnt
    %% first tactile/auditory stimuli
    lfirstAudStimR{iTrials} = false(frames, fsPostTime);
    rfirstAudStimR{iTrials} = false(frames, fsPostTime);
    lfirstTacStimR{iTrials} = false(frames, fsPostTime);
    rfirstTacStimR{iTrials} = false(frames, fsPostTime);

    firstStim = NaN;
    if bhv.StimType(iTrials) == 2 || bhv.StimType(iTrials) == 6 %auditory or mixed stimulus
        if ~isempty(audStimL{iTrials}(~isnan(audStimL{iTrials})))
            firstStim = round((audStimL{iTrials}(1) + preStimDur) * sRate);
            stimEnd = firstStim - 1 + fsPostTime; stimEnd = min([frames stimEnd]);
            lfirstAudStimR{iTrials}(:,1 : stimEnd - firstStim + 1)  = timeR(:, firstStim : stimEnd);
        end
        if ~isempty(audStimR{iTrials}(~isnan(audStimR{iTrials})))
            firstStim = round((audStimR{iTrials}(1) + preStimDur) * sRate);
            stimEnd = firstStim - 1 + fsPostTime; stimEnd = min([frames stimEnd]);
            rfirstAudStimR{iTrials}(:,1 : stimEnd - firstStim + 1) = timeR(:, firstStim : stimEnd);
        end
    end

    if bhv.StimType(iTrials) == 4 || bhv.StimType(iTrials) == 6 %tactile or mixed stimulus
        if ~isempty(tacStimL{iTrials}(~isnan(tacStimL{iTrials})))
            firstStim = round((tacStimL{iTrials}(1) + preStimDur) * sRate);
            stimEnd = firstStim - 1 + fsPostTime; stimEnd = min([frames stimEnd]);
            lfirstTacStimR{iTrials}(:,1 : stimEnd - firstStim + 1) = timeR(:, firstStim : stimEnd);
        end
        if ~isempty(tacStimR{iTrials}(~isnan(tacStimR{iTrials})))
            firstStim = round((tacStimR{iTrials}(1) + preStimDur) * sRate);
            stimEnd = firstStim - 1 + fsPostTime; stimEnd = min([frames stimEnd]);
            rfirstTacStimR{iTrials}(:,1 : stimEnd - firstStim + 1) = timeR(:, firstStim : stimEnd);
        end
    end

    if gaussShift > 1
        % subsample regressors
        lfirstTacStimR{iTrials} = lfirstTacStimR{iTrials}(:,1:gaussShift:end);
        rfirstTacStimR{iTrials} = rfirstTacStimR{iTrials}(:,1:gaussShift:end);
        lfirstAudStimR{iTrials} = lfirstAudStimR{iTrials}(:,1:gaussShift:end);
        rfirstAudStimR{iTrials} = rfirstAudStimR{iTrials}(:,1:gaussShift:end);
    end

    %% other tactile/auditory stimuli
    lAudStimR{iTrials} = false(frames, sPostTime);
    rAudStimR{iTrials} = false(frames, sPostTime);
    lTacStimR{iTrials} = false(frames, sPostTime);
    rTacStimR{iTrials} = false(frames, sPostTime);

    for iRegs = 0 : sPostTime - 1
        allStims = audStimL{iTrials}(2:end) + (iRegs * 1/sRate);
        lAudStimR{iTrials}(logical(histcounts(allStims,-preStimDur:1/sRate:postStimDur)),iRegs+1) = 1;

        allStims = audStimR{iTrials}(2:end) + (iRegs * 1/sRate);
        rAudStimR{iTrials}(logical(histcounts(allStims,-preStimDur:1/sRate:postStimDur)),iRegs+1) = 1;

        allStims = tacStimL{iTrials}(2:end) + (iRegs * 1/sRate);
        lTacStimR{iTrials}(logical(histcounts(allStims,-preStimDur:1/sRate:postStimDur)),iRegs+1) = 1;

        allStims = tacStimR{iTrials}(2:end) + (iRegs * 1/sRate);
        rTacStimR{iTrials}(logical(histcounts(allStims,-preStimDur:1/sRate:postStimDur)),iRegs+1) = 1;
    end

    if gaussShift > 1
        % subsample regressors
        lTacStimR{iTrials} = lTacStimR{iTrials}(:,1:gaussShift:end);
        rTacStimR{iTrials} = rTacStimR{iTrials}(:,1:gaussShift:end);
        lAudStimR{iTrials} = lAudStimR{iTrials}(:,1:gaussShift:end);
        rAudStimR{iTrials} = rAudStimR{iTrials}(:,1:gaussShift:end);
    end

    %% spout regressors
    spoutIdx = round((preStimDur + spoutTime(iTrials)) * sRate) : round((preStimDur + postStimDur) * sRate); %index for which part of the trial should be covered by spout regressors
    spoutR{iTrials} = false(frames, maxSpoutRegs);
    if ~isnan(spoutTime(iTrials))
        spoutR{iTrials}(:, 1:length(spoutIdx)) = timeR(:, spoutIdx);
    end

    spoutOutR{iTrials} = false(frames, 3);
    spoutOut = round((preStimDur + spoutOutTime(iTrials)) * sRate); %time when opposing spout moved out again
    if ~isnan(spoutOut) && spoutOut < (frames + 1)
        cInd = spoutOut : spoutOut + 2; cInd(cInd > frames) = [];
        temp = diag(ones(1,3));
        spoutOutR{iTrials}(cInd, :) = temp(1:length(cInd),:);
    end

    if gaussShift > 1
        % subsample regressors
        spoutR{iTrials} = spoutR{iTrials}(:,1:gaussShift:end);
        spoutOutR{iTrials} = spoutOutR{iTrials}(:,1:gaussShift:end);
    end

    %% lick regressors
    lLickR{iTrials} = false(frames, length(motorIdx));
    rLickR{iTrials} = false(frames, length(motorIdx));

    for iRegs = 0 : length(motorIdx)-1
        licks = lickL{iTrials} - ((mPreTime/sRate) - (iRegs * 1/sRate));
        lLickR{iTrials}(logical(histcounts(licks,-preStimDur:1/sRate:postStimDur)),iRegs+1) = 1;

        licks = lickR{iTrials} - ((mPreTime/sRate) - (iRegs * 1/sRate));
        rLickR{iTrials}(logical(histcounts(licks,-preStimDur:1/sRate:postStimDur)),iRegs+1) = 1;
    end

    if gaussShift > 1
        % subsample regressors
        lLickR{iTrials} = lLickR{iTrials}(:,1:gaussShift:end);
        rLickR{iTrials} = rLickR{iTrials}(:,1:gaussShift:end);
    end

    %% lever in
    leverInR{iTrials} = false(frames, leverMoveDur);
    leverShift = round((preStimDur + leverIn(iTrials))* sRate); %timepoint in frames when lever moved in, relative to lever grab

    if ~isnan(leverShift)
        if leverShift > 0 %lever moved in during the recorded trial
            leverInR{iTrials}(leverShift : leverShift + leverMoveDur -1, :) = diag(ones(1, leverMoveDur));
        elseif (leverShift + leverMoveDur) > 0  %lever was moving before data was recorded but still moving at trial onset
            leverInR{iTrials}(1 : leverMoveDur + leverShift, :) = [zeros(leverMoveDur + leverShift, abs(leverShift)) diag(ones(1, leverMoveDur + leverShift))];
        end
    end

    if gaussShift > 1
        % subsample regressors
        leverInR{iTrials} = leverInR{iTrials}(:,1:gaussShift:end);
    end

    %% dual-handle indicator sound
    handleSoundR{iTrials} = false(frames, sPostTime);
    for iRegs = 0 : sPostTime - 1
        allStims = handleSounds{iTrials}(1:end) + (iRegs * 1/sRate);
        allStims = allStims(~isnan(allStims)) + 0.001;
        handleSoundR{iTrials}(logical(histcounts(allStims,-preStimDur:1/sRate:postStimDur)),iRegs+1) = 1;
    end

    %% choice, reward, and state
    stimTemp = false(frames, frames + maxStimShift);
    stimShift = round(stimTime(iTrials) * sRate); %amount of stimshift compared to possible maximum. move diagonal on x-axis accordingly.

    if (stimShift > maxStimShift) || isnan(stimTime(iTrials))
        stimTemp = NaN(frames, frames + maxStimShift); %don't use trial if stim onset is too late
    else
        stimTemp(:, end - stimShift - frames + 1 : end - stimShift) = timeR;
    end

    rewardR{iTrials} = false(size(stimTemp));

    if bhv.Rewarded(iTrials) %rewarded
        rewardR{iTrials} = stimTemp; %trial was rewarded
    end

    % state as binary design matrix
    stateR{iTrials} = false(size(stimTemp));
    if states(iTrials) %if engaged
        stateR{iTrials} = stimTemp;
    end

    % get L/R choices as binary design matrix
    ChoiceR{iTrials} = false(size(stimTemp));
    if bhv.ResponseSide(iTrials) == 1 %IF LEFT CHOICE!!! Left choice is 1, right is 2
        ChoiceR{iTrials} = stimTemp;
    end


    % previous trial regressors
    if iTrials == 1 %don't use first trial
        prevRewardR{iTrials} = NaN(size(timeR(:,1:end-4)));
        prevChoiceR{iTrials} = NaN(size(timeR(:,1:end-4)));
        prevStimR{iTrials} = NaN(size(timeR(:,1:end-4)));

    else %for all subsequent trials
        % same as for regular choice regressors but for prevoious trial
        prevChoiceR{iTrials} = false(size(timeR(:,1:end-4)));
        if SessionData.ResponseSide(bTrials(iTrials)-1) == 1
            prevChoiceR{iTrials} = timeR(:,1:end-4);
        end

        prevStimR{iTrials} = false(size(timeR(:,1:end-4)));
        if SessionData.CorrectSide(bTrials(iTrials)-1) == 1 % if previous trial had a left target
            prevStimR{iTrials} = timeR(:,1:end-4);
        end

        prevRewardR{iTrials} = false(size(timeR(:,1:end-4)));
        if SessionData.Rewarded(bTrials(iTrials)-1) %last trial was rewarded
            prevRewardR{iTrials} = timeR(:,1:end-4);
        end
    end

    % subsequent trial regressors
    if iTrials == length(bTrials) %don't use lat trial
        nextChoiceR{iTrials} = NaN(size(timeR(:,1:end-4)));
        repeatChoiceR{iTrials} = NaN(size(timeR(:,1:end-4)));

    else %for all subsequent trials
        nextChoiceR{iTrials} = false(size(timeR(:,1:end-4)));
        if SessionData.ResponseSide(bTrials(iTrials)+1) == 1 %choice in next trial is left
            nextChoiceR{iTrials} = timeR(:,1:end-4);
        end

        repeatChoiceR{iTrials} = false(size(timeR(:,1:end-4)));
        if SessionData.ResponseSide(bTrials(iTrials)) == SessionData.ResponseSide(bTrials(iTrials)+1) %choice in next trial is similar to the current one
            repeatChoiceR{iTrials} = timeR(:,1:end-4);
        end
    end

    if gaussShift > 1
        % subsample regressors
        rewardR{iTrials} = rewardR{iTrials}(:,1:gaussShift:end);
        prevRewardR{iTrials} = prevRewardR{iTrials}(:,1:gaussShift:end);
        ChoiceR{iTrials} = ChoiceR{iTrials}(:,1:gaussShift:end);
        prevChoiceR{iTrials} = prevChoiceR{iTrials}(:,1:gaussShift:end);
        prevStimR{iTrials} = prevStimR{iTrials}(:,1:Shift:end);
        nextChoiceR{iTrials} = nextChoiceR{iTrials}(:,1:gaussShift:end);
        repeatChoiceR{iTrials} = repeatChoiceR{iTrials}(:,1:gaussShift:end);
    end

    %determine timepoint of reward given
    waterR{iTrials} = false(frames, sRate * 2);
    if ~isnan(water(iTrials)) && ~isempty(water(iTrials))
        waterOn = round((preStimDur + water(iTrials)) * sRate); %timepoint in frames when reward was given
        if waterOn <= frames
            waterR{iTrials}(:, 1 : size(timeR,2) - waterOn + 1) = timeR(:, waterOn:end);
        end
    end

    if gaussShift > 1
        waterR{iTrials} = waterR{iTrials}(:,1:gaussShift:end); % subsample regressor
    end

    %% lever grabs
    cGrabs = levGrabL{iTrials};
    cGrabs(cGrabs >= postStimDur) = []; %remove grabs after end of imaging
    cGrabs(find(diff(cGrabs) < tapDur) + 1) = []; %remove grabs that are too close to one another
    lGrabR{iTrials} = histcounts(cGrabs,-preStimDur:1/sRate:postStimDur)'; %convert to binary trace

    cGrabs = levGrabR{iTrials};
    cGrabs(cGrabs >= postStimDur) = []; %remove grabs after end of imaging
    cGrabs(find(diff(cGrabs) < tapDur) + 1) = []; %remove grabs that are too close to one another
    rGrabR{iTrials} = histcounts(cGrabs,-preStimDur:1/sRate:postStimDur)'; %convert to binary trace

    %% pupil / whisk / nose / face / body regressors
    bhvFrameRate = round(1/mean(diff(pTime{bTrials(iTrials)}))); %framerate of face camera
    trialOn = bhv.TrialStartTime(iTrials) + (stimGrab(iTrials) - preStimDur);
    trialTime = pTime{bTrials(iTrials)} - trialOn;
    rejIdx = trialTime < trialDur; %don't use late frames
    trialTime = trialTime(rejIdx);

    if isempty(trialTime) || trialTime(1) > 0 %check if there is missing time at the beginning of a trial
        warning(['Trial ' int2str(bTrials(iTrials)) ': Missing behavioral video frames at trial onset. Trial removed from analysis']);
        fastPupilR{iTrials} = NaN(frames, 1);
        slowPupilR{iTrials} = NaN(frames, 1);
        whiskR{iTrials} = NaN(frames, 1);
        noseR{iTrials} = NaN(frames, 1);
        faceR{iTrials} = NaN(frames, 1);
        bodyR{iTrials} = NaN(frames, 1);

    else
        timeLeft = trialDur - trialTime(end); %check if there is missing time at the end of a trial
        if (timeLeft < trialDur * 0.9) && (timeLeft > 0) %if there is some time missing to make a whole trial
            addTime = trialTime(end) + (1/bhvFrameRate : 1/bhvFrameRate : timeLeft + 1/bhvFrameRate); %add some dummy times to make complete trial
            trialTime = [trialTime' addTime];
        end

        fastPupilR{iTrials} = Behavior_vidResamp(fPupil{bTrials(iTrials)}(rejIdx), trialTime, sRate);
        fastPupilR{iTrials} = smooth(fastPupilR{iTrials}(end - frames + 1 : end), 'rlowess');

        slowPupilR{iTrials} = Behavior_vidResamp(sPupil{bTrials(iTrials)}(rejIdx), trialTime, sRate);
        slowPupilR{iTrials} =  smooth(slowPupilR{iTrials}(end - frames + 1 : end), 'rlowess');

        whiskR{iTrials} = Behavior_vidResamp(whisker{bTrials(iTrials)}(rejIdx), trialTime, sRate);
        whiskR{iTrials} = smooth(whiskR{iTrials}(end - frames + 1 : end), 'rlowess');

        noseR{iTrials} = Behavior_vidResamp(nose{bTrials(iTrials)}(rejIdx), trialTime, sRate);
        noseR{iTrials} = smooth(noseR{iTrials}(end - frames + 1 : end), 'rlowess');

        faceR{iTrials} = Behavior_vidResamp(faceM{bTrials(iTrials)}(rejIdx), trialTime, sRate);
        faceR{iTrials} = smooth(faceR{iTrials}(end - frames + 1 : end), 'rlowess');


        % body regressors
        bhvFrameRate = round(1/mean(diff(bTime{bTrials(iTrials)}))); %framerate of body camera
        trialTime = bTime{bTrials(iTrials)} - trialOn;
        rejIdx = trialTime < trialDur; %don't use late frames
        trialTime = trialTime(rejIdx);
        timeLeft = trialDur - trialTime(end); %check if there is missing time at the end of a trial

        if (timeLeft < trialDur * 0.9) && (timeLeft > 0) %if there is some time missing to make a whole trial
            addTime = trialTime(end) + (1/bhvFrameRate : 1/bhvFrameRate : timeLeft + 1/bhvFrameRate); %add some dummy times to make complete trial
            trialTime = [trialTime' addTime];
        end

        bodyR{iTrials} = Behavior_vidResamp(bodyM{bTrials(iTrials)}(rejIdx), trialTime, sRate);
        bodyR{iTrials} = smooth(bodyR{iTrials}(end - frames + 1 : end), 'rlowess');
    end

    %% piezo sensor information
    if strcmpi(dType,'Widefield')
        if exist([cPath 'Analog_'  num2str(trials(iTrials)) '.dat'],'file') ~= 2  %check if files exists on hdd and pull from server otherwise
            cFile = dir([sPath 'Analog_'  num2str(trials(iTrials)) '.dat']);
            copyfile([sPath 'Analog_'  num2str(trials(iTrials)) '.dat'],[cPath 'Analog_'  num2str(trials(iTrials)) '.dat']);
        end
        [~,Analog] = Widefield_LoadData([cPath 'Analog_'  num2str(trials(iTrials)) '.dat'],'Analog'); %load analog data
        stimOn = find(diff(double(Analog(stimLine,:)) > 1500) == 1); %find stimulus onset in current trial
    elseif strcmpi(dType,'twoP')
        Analog = squeeze(data.analog(:,:,iTrials));
        stimOn = find(diff(double(Analog(stimLine,:)) > 1) == 1); %find stimulus onset in current trial
    end

    if ~isnan(stimTime(iTrials))
        try
            Analog(1,round(stimOn + ((postStimDur-stimTime(iTrials)) * 1000) - 1)) = 0; %make sure there are enough datapoints in analog signal
            temp = Analog(piezoLine,round(stimOn - ((preStimDur + stimTime(iTrials)) * 1000)) : round(stimOn + ((postStimDur - stimTime(iTrials))* 1000) - 1)); % data from piezo sensor. Should encode animals hindlimb motion.
            temp = smooth(double(temp), sRate*5, 'lowess')'; %do some smoothing
            temp = [repmat(temp(1),1,1000) temp repmat(temp(end),1,1000)]; %add some padding on both sides to avoid edge effects when resampling
            temp = resample(double(temp), sRate, 1000); %resample to imaging rate
            piezoR{iTrials} = temp(sRate + 1 : end - sRate)'; %remove padds again
            piezoR{iTrials} = piezoR{iTrials}(end - frames + 1:end); %make sure, the length is correct

            temp = abs(hilbert(diff(piezoR{iTrials})));
            piezoMoveR{iTrials} = [temp(1); temp]; %keep differential motion signal
            clear temp
        catch
            piezoMoveR{iTrials} = NaN(frames, 1);
            piezoR{iTrials} = NaN(frames, 1);
        end
    else
        piezoMoveR{iTrials} = NaN(frames, 1);
        piezoR{iTrials} = NaN(frames, 1);
    end

    % give some feedback over progress
    if rem(iTrials,50) == 0
        fprintf(1, 'Current trial is %d out of %d\n', iTrials,trialCnt);
        toc
    end
end

%% get proper design matrices for handle grab
lGrabR = cat(1,lGrabR{:});
lGrabR = Widefield_analogToDesign(lGrabR, 0.5, trialCnt, sRate, sRate, motorIdx, gaussShift); %get design matrix

rGrabR = cat(1,rGrabR{:});
rGrabR = Widefield_analogToDesign(rGrabR, 0.5, trialCnt, sRate, sRate, motorIdx, gaussShift); %get design matrix

%% rebuild analog motor regressors to get proper design matrices
temp = double(cat(1,fastPupilR{:}));
temp = (temp - prctile(temp,1))./ nanstd(temp); %minimum values are at 0, signal in standard deviation units
[dMat, traceOut] = Widefield_analogToDesign(temp, median(temp), trialCnt, sRate, sRate, motorIdx, gaussShift);
[dMat2, ~] = Widefield_analogToDesign(traceOut, prctile(traceOut,75), trialCnt, sRate, sRate, motorIdx, gaussShift);

fastPupilAnalog = single(traceOut);
fastPupilDigital = cat(1,dMat2{:});
%fastPupilHiDigital = cat(1,dMat2{:});

[dMat, traceOut] = Widefield_analogToDesign(double(cat(1,whiskR{:})), 0.5, trialCnt, sRate, sRate, motorIdx, gaussShift);
[dMat2, ~] = Widefield_analogToDesign(double(cat(1,whiskR{:})), 2, trialCnt, sRate, sRate, motorIdx, gaussShift);

whiskAnalog = single(traceOut);
whiskDigital = cat(1,dMat2{:});
%whiskHiDigital = cat(1,dMat2{:});

[dMat, traceOut] = Widefield_analogToDesign(double(cat(1,noseR{:})), 0.5, trialCnt, sRate, sRate, motorIdx, gaussShift);
[dMat2, ~] = Widefield_analogToDesign(double(cat(1,noseR{:})), 2, trialCnt, sRate, sRate, motorIdx, gaussShift);

noseAnalog = single(traceOut);
noseDigital = cat(1,dMat2{:});
%noseHiDigital = cat(1,dMat2{:});

[dMat, traceOut] = Widefield_analogToDesign(double(cat(1,piezoR{:})), 2, trialCnt, sRate, sRate, motorIdx, gaussShift);
[dMat2, traceOut2] = Widefield_analogToDesign(double(cat(1,piezoMoveR{:})), 0.5, trialCnt, sRate, sRate, motorIdx, gaussShift);
[dMat3, ~] = Widefield_analogToDesign(double(cat(1,piezoMoveR{:})), 2, trialCnt, sRate, sRate, motorIdx, gaussShift);

piezoAnalog = traceOut;
piezoDigital = cat(1,dMat{:});
piezoMoveAnalog = traceOut2;
piezoMoveDigital = cat(1,dMat3{:});
%piezoMoveHiDigital = cat(1,dMat3{:});

[dMat, traceOut] = Widefield_analogToDesign(double(cat(1,faceR{:})), 0.5, trialCnt, sRate, sRate, motorIdx, gaussShift);
[dMat2, ~] = Widefield_analogToDesign(double(cat(1,faceR{:})), 2, trialCnt, sRate, sRate, motorIdx, gaussShift);

faceAnalog = single(traceOut);
faceDigital = cat(1,dMat2{:});
%faceHiDigital = cat(1,dMat2{:});

[dMat, traceOut] = Widefield_analogToDesign(double(cat(1,bodyR{:})), 0.5, trialCnt, sRate, sRate, motorIdx, gaussShift);
[dMat2, ~] = Widefield_analogToDesign(double(cat(1,bodyR{:})), 2, trialCnt, sRate, sRate, motorIdx, gaussShift);

bodyAnalog = single(traceOut);
bodyDigital = cat(1,dMat2{:});
%bodyHiDigital = cat(1,dMat2{:});

clear piezoR1 piezoR2 dMat dMat2 traceOut temp

%% re-align behavioral video data and Vc to lever grab instead of stimulus onset
if strcmpi(dType,'Widefield')
    iiSpikeFrames = findInterictalSpikes(U, Vc, 2, false); %find interictal spikes
    Vc = interpOverInterictal(Vc, iiSpikeFrames); %interpolate over interictal spikes
end

V1 = reshape(V1, [], SessionData.nTrials, bhvDimCnt);
V2 = reshape(V2, [], SessionData.nTrials, bhvDimCnt);

%if video sampling rate is different from widefield, resample video data
if bhvOpts.targRate ~= sRate
    vidR = NaN(size(Vc,2), length(bTrials), size(V1,3), 'single');
    moveR = NaN(size(Vc,2), length(bTrials), size(V1,3), 'single');
    for iTrials = 1 : length(bTrials)

        temp1 = squeeze(V1(:,bTrials(iTrials),:));
        if ~any(isnan(temp1(:)))
            trialTime = 1/bhvRate : 1/bhvRate : size(Vc,2)/sRate;
            vidR(:, iTrials, :) = Behavior_vidResamp(double(temp1), trialTime, sRate);

            temp2 = squeeze(V2(:,bTrials(iTrials),:));
            trialTime = 1/bhvRate : 1/bhvRate : size(Vc,2)/sRate;
            moveR(:, iTrials, :) = Behavior_vidResamp(double(temp2), trialTime, sRate);
        end
    end
else
    vidR = V1(:,bTrials,:); clear V1 %get correct trials from behavioral video data.
    moveR = V2(:,bTrials,:); clear V2 %get correct trials from behavioral video data.
end

% re-align video data
temp1 = NaN(dims,frames,trialCnt);
temp2 = NaN(frames,trialCnt,bhvDimCnt);
temp3 = NaN(frames,trialCnt,bhvDimCnt);
temp4 = NaN(2,frames,trialCnt);
for x = 1 : size(vidR,2) %iterate thru trials
    try
        temp1(:,:,x) = Vc(:,(shVal - ceil(stimTime(x) / (1/sRate))) - (preStimDur * sRate) : (shVal - ceil(stimTime(x) / (1/sRate))) + (postStimDur * sRate) - 1,x);
        temp2(:,x,:) = vidR((shVal - ceil(stimTime(x) / (1/sRate))) - (preStimDur * sRate) : (shVal - ceil(stimTime(x) / (1/sRate))) + (postStimDur * sRate) - 1,x,:);
        temp3(:,x,:) = moveR((shVal - ceil(stimTime(x) / (1/sRate))) - (preStimDur * sRate) : (shVal - ceil(stimTime(x) / (1/sRate))) + (postStimDur * sRate) - 1,x,:);
        if strcmpi(dType,'twoP')
            temp4(:,:,x) = data.DS(:,(shVal - ceil(stimTime(x) / (1/sRate))) - (preStimDur * sRate) : (shVal - ceil(stimTime(x) / (1/sRate))) + (postStimDur * sRate) - 1,x);
        end
    catch
        fprintf(1,'Could not align trial %d. Relative stim time: %fs\n', x, stimTime(x));
    end
end
Vc = reshape(temp1,dims,[]); clear temp1
vidR = reshape(temp2,[],bhvDimCnt); clear temp2
moveR = reshape(temp3,[],bhvDimCnt); clear temp3

if strcmpi(dType,'twoP')
    DS = reshape(temp4,2,[]); %keep image motion trace for 2p imaging
end
clear temp4


%% reshape regressors, make design matrix and indices for regressors that are used for the model
timeR = repmat(logical(diag(ones(1,frames))),trialCnt,1); %time regressor
timeR = timeR(:,1:end-4);

lGrabR = cat(1,lGrabR{:});
lGrabRelR = cat(1,lGrabRelR{:}); %not in the model currently
rGrabR = cat(1,rGrabR{:});
rGrabRelR = cat(1,rGrabRelR{:}); %not in the model currently

lLickR = cat(1,lLickR{:});
rLickR = cat(1,rLickR{:});
leverInR = cat(1,leverInR{:});
leverInR(:,sum(leverInR) == 0) = [];

handleSoundR = cat(1,handleSoundR{:});

lTacStimR = cat(1,lTacStimR{:});
rTacStimR = cat(1,rTacStimR{:});
lAudStimR = cat(1,lAudStimR{:});
rAudStimR = cat(1,rAudStimR{:});

lfirstTacStimR = cat(1,lfirstTacStimR{:});
rfirstTacStimR = cat(1,rfirstTacStimR{:});
lfirstAudStimR = cat(1,lfirstAudStimR{:});
rfirstAudStimR = cat(1,rfirstAudStimR{:});

spoutR = cat(1,spoutR{:});
spoutOutR = cat(1,spoutOutR{:});
spoutR(:,sum(spoutR) == 0) = [];
spoutOutR(:,sum(spoutOutR) == 0) = [];

rewardR = cat(1,rewardR{:});
prevRewardR = cat(1,prevRewardR{:});

ChoiceR = cat(1,ChoiceR{:});

stateR = cat(1,stateR{:});

prevChoiceR = cat(1,prevChoiceR{:});
prevStimR = cat(1,prevStimR{:});
nextChoiceR = cat(1,nextChoiceR{:});
repeatChoiceR = cat(1,repeatChoiceR{:});

waterR = cat(1,waterR{:});

slowPupilR = cat(1,slowPupilR{:});
slowPupilR(~isnan(slowPupilR(:,1)),:) = zscore(slowPupilR(~isnan(slowPupilR(:,1)),:));

%% create full design matrix

fullR = {timeR ChoiceR rewardR lGrabR rGrabR  ...
    lLickR rLickR handleSoundR lfirstTacStimR lTacStimR rfirstTacStimR rTacStimR ...
    lfirstAudStimR lAudStimR rfirstAudStimR rAudStimR prevRewardR prevChoiceR ...
    nextChoiceR waterR piezoAnalog piezoDigital piezoMoveAnalog piezoMoveDigital  whiskAnalog whiskDigital  noseAnalog noseDigital  fastPupilAnalog fastPupilDigital  ...
    slowPupilR faceAnalog faceDigital  bodyAnalog bodyDigital  moveR vidR stateR};



% labels for different regressor sets. It is REALLY important this agrees with the order of regressors in fullR.
regLabels = {
    'time' 'Choice' 'reward' 'lGrab' 'rGrab' 'lLick' 'rLick' 'handleSound' ...
    'lfirstTacStim' 'lTacStim' 'rfirstTacStim' 'rTacStim' 'lfirstAudStim' 'lAudStim' 'rfirstAudStim' 'rAudStim' ...
    'prevReward' 'prevChoice' 'nextChoice' 'water' 'piezoAnalog' 'piezoDigital' 'piezoMoveAnalog' 'piezoMoveDigital' ...
    'whiskAnalog' 'whiskDigital' 'noseAnalog' 'noseDigital' 'fastPupilAnalog' 'fastPupilDigital' 'slowPupil' 'faceAnalog' ...
    'faceDigital' 'bodyAnalog' 'bodyDigital' 'Move' 'bhvVideo' 'state'};

%index to reconstruct different response kernels
regIdx = [];
for i = 1:length(regLabels)
    t = ones(1,size(fullR{i},2))*find(ismember(regLabels,regLabels{i}));
    regIdx = [regIdx, t];
end

fullR = cat(2,fullR{:});
%% compute the indices of the event that we're aligning to - for easier recovery of kernels

preStimFrames = preStimDur*sRate;
postStimFrames = postStimDur*sRate;

timeZero = zeros(1,size(timeR,2)); timeZero(preStimFrames) = 1; %handle grab time
choiceZero = zeros(1,size(ChoiceR,2)); choiceZero(preStimFrames + maxStimShift) = 1; %stim on time choice is realigned to stimulus
stateZero = zeros(1,size(stateR,2)); stateZero(preStimFrames + maxStimShift) = 1;
rewardZero = zeros(1,size(rewardR,2)); rewardZero(preStimFrames + maxStimShift) = 1; %stim on time, same alignment as choice
lGrabZero = zeros(1,size(lGrabR,2)); lGrabZero(mPreTime + 1) = 1; %frame of actual motor event
rGrabZero = zeros(1,size(rGrabR,2)); rGrabZero(mPreTime + 1) = 1;
lLickZero = zeros(1,size(lLickR,2)); lLickZero(mPreTime + 1) = 1;
rLickZero = zeros(1,size(rLickR,2)); rLickZero(mPreTime + 1) = 1;
handleSoundZero = zeros(1,size(handleSoundR,2)); handleSoundZero(1) = 1; %first frame marks event
lfirstTacStimZero = zeros(1,size(lfirstTacStimR,2)); lfirstTacStimZero(1) = 1; %first frame marks event
lTacStimZero = zeros(1,size(lTacStimR,2)); lTacStimZero(1) = 1;
rfirstTacStimZero = zeros(1,size(rfirstTacStimR,2)); rfirstTacStimZero(1) = 1;
rTacStimZero = zeros(1,size(rTacStimR,2)); rTacStimZero(1) = 1;
lfirstAudStimZero = zeros(1,size(lfirstAudStimR,2)); lfirstAudStimZero(1) = 1; %first frame marks event
lAudStimZero = zeros(1,size(lAudStimR,2)); lAudStimZero(1) = 1;
rfirstAudStimZero = zeros(1,size(rfirstAudStimR,2)); rfirstAudStimZero(1) = 1;
rAudStimZero = zeros(1,size(rAudStimR,2)); rAudStimZero(1) = 1;
prevRewardZero = zeros(1,size(prevRewardR,2)); prevRewardZero(preStimFrames) = 1; %previous trial regressors are aligned to handle grab on current trial
prevChoiceZero = zeros(1,size(prevChoiceR,2)); prevChoiceZero(preStimFrames) = 1;
nextChoiceZero = zeros(1,size(nextChoiceR,2)); nextChoiceZero(preStimFrames) = 1;
waterZero = zeros(1,size(waterR,2)); waterZero(1) = 1; %begins at water delivery, lots of trials are cut off by trial end
piezoAnalogZero = zeros(1,size(piezoAnalog,2)); %analog regressors have no alignment value
piezoDigitalZero = zeros(1,size(piezoDigital,2)); piezoDigitalZero(mPreTime + 1) = 1;
piezoMoveAnalogZero = zeros(1,size(piezoMoveAnalog,2));
piezoMoveDigitalZero = zeros(1,size(piezoMoveDigital,2)); piezoMoveDigitalZero(mPreTime + 1) = 1;
%piezoMoveHiDigitalZero = zeros(1,size(piezoMoveHiDigital,2)); piezoMoveHiDigitalZero(mPreTime + 1) = 1;
whiskAnalogZero = zeros(1,size(whiskAnalog,2));
whiskDigitalZero = zeros(1,size(whiskDigital,2)); whiskDigitalZero(mPreTime + 1) = 1;
%whiskHiDigitalZero = zeros(1,size(whiskHiDigital,2)); whiskHiDigitalZero(mPreTime + 1) = 1;
noseAnalogZero = zeros(1,size(noseAnalog,2));
noseDigitalZero = zeros(1,size(noseDigital,2)); noseDigitalZero(mPreTime + 1) = 1;
%noseHiDigitalZero = zeros(1,size(noseHiDigital,2)); noseHiDigitalZero(mPreTime + 1) = 1;
fastPupilAnalogZero = zeros(1,size(fastPupilAnalog,2));
fastPupilDigitalZero = zeros(1,size(fastPupilDigital,2)); fastPupilDigitalZero(mPreTime + 1) = 1;
%fastPupilHiDigitalZero = zeros(1,size(fastPupilHiDigital,2)); fastPupilHiDigitalZero(mPreTime + 1) = 1;
slowPupilZero = zeros(1,size(slowPupilR,2));
faceAnalogZero = zeros(1,size(faceAnalog,2));
faceDigitalZero = zeros(1,size(faceDigital,2)); faceDigitalZero(mPreTime + 1) = 1;
%faceHiDigitalZero = zeros(1,size(faceHiDigital,2)); faceHiDigitalZero(mPreTime + 1) = 1;
bodyAnalogZero = zeros(1,size(bodyAnalog,2));
bodyDigitalZero = zeros(1,size(bodyDigital,2)); bodyDigitalZero(mPreTime + 1) = 1;
%bodyHiDigitalZero = zeros(1,size(bodyHiDigital,2)); bodyHiDigitalZero(mPreTime + 1) = 1;
moveRZero = zeros(1,size(moveR,2));
vidRZero = zeros(1,size(vidR,2));

regZeroFrames = [timeZero choiceZero rewardZero lGrabZero rGrabZero lLickZero rLickZero handleSoundZero ...
    lfirstTacStimZero lTacStimZero rfirstTacStimZero rTacStimZero lfirstAudStimZero lAudStimZero ...
    rfirstAudStimZero rAudStimZero prevRewardZero prevChoiceZero nextChoiceZero waterZero piezoAnalogZero ...
    piezoDigitalZero piezoMoveAnalogZero  piezoMoveDigitalZero whiskAnalogZero ...
    whiskDigitalZero noseAnalogZero noseDigitalZero fastPupilAnalogZero ...
    fastPupilDigitalZero slowPupilZero faceAnalogZero faceDigitalZero ...
    bodyAnalogZero bodyDigitalZero moveRZero vidRZero stateZero];

%%

% orthogonalize video against spout/handle movement
vidIdx = find(ismember(regIdx, find(ismember(regLabels,{'Move' 'bhvVideo'})))); %index for video regressors
trialIdx = ~isnan(mean(fullR(:,vidIdx),2)); %don't use trials that failed to contain behavioral video data

if orthogvideotoinstructed
    smallR = [lGrabR rGrabR lLickR rLickR leverInR spoutR spoutOutR];
else
    smallR = [leverInR spoutR spoutOutR];
end
%temp = fullR;
for iRegs = 1 : length(vidIdx)
    matrix_to_orthogonalize = [smallR(trialIdx,:) fullR(trialIdx,vidIdx(iRegs))];
    [Q,~] = qr(matrix_to_orthogonalize,0); %orthogonalize video against other regressors
    %Q = qr(matrix_to_orthogonalize,0); %orthogonalize video against other regressors
    fullR(trialIdx,vidIdx(iRegs)) = Q(:,end); % transfer orthogonolized video regressors back to design matrix
end
%figure;imshow(temp(1:2000,:))
%figure;imshow(fullR(1:2000,:))


% reject trials with broken regressors that contain NaNs
trialIdx = isnan(mean(fullR,2)); %don't use first trial or trials that failed to contain behavioral video data
fprintf(1, 'Rejected %d/%d trials for NaN entries in regressors\n', sum(trialIdx)/frames, trialCnt);
fullR(trialIdx,:) = []; %clear bad trials

Vc(:,trialIdx) = []; %clear bad trials
zeromeanVc = bsxfun(@minus, Vc, mean(Vc,2)); %should be zero-mean

trialIdx = trialIdx(1:75:end);
usedTrials = usedTrials(~trialIdx);
assert(length(usedTrials) == size(fullR,1) / frames, 'there is a mismatch between usedTrials and the length of the design matrix')

end


