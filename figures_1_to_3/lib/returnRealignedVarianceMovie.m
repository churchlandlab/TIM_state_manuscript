function alignedVariance = returnRealignedVarianceMovie(cpath,mouse,rec,modelfile,segIdx,nframes)

PRESTIMDUR = 2;
POSTSTIMDUR = 3;
SRATE = 15;
HANDLEGRABFRAME = 30;
fprintf('Realignment assumes that the handle grab is on frame %i\n',HANDLEGRABFRAME);



segFrames = cumsum(floor(segIdx * SRATE));
%%
datapath = [cpath filesep mouse filesep 'SpatialDisc' filesep rec filesep];
load('C:\Data\churchland\ridgeModel\widefield\allenDorsalMapSM.mat','dorsalMaps'); %Get allen atlas
try
    load([datapath modelfile],'zeromeanVc','Vm', 'U','usedTrials');
catch
    fprintf('\nThe encoding model file does not exist!');
    alignedVariance = NaN;
    return
end

%encoding model data is currently aligned to handle grab, we need to
%realign it to different task epochs for better visualization
f = dir([datapath filesep '*SpatialDisc*.mat']);
SessionData = []; %needed to avoid static workspace error
load([f.folder filesep f.name]); clear f;
cBhv = selectBehaviorTrials(SessionData,usedTrials); %select trials that were used for the encoding model

zeromeanVc = reshape(zeromeanVc,200,nframes,[]);
Vm = reshape(Vm,200,nframes,[]);

alignedVc = innerAlignFunc(zeromeanVc);
alignedVm = innerAlignFunc(Vm);

alingnednframes = size(alignedVc,2);
alignedVc = reshape(alignedVc,200,[]);
alignedVm = reshape(alignedVm,200,[]);

% movie for predicted variance - after realignment
alignedVariance = computeFramewiseVariance(alignedVc, alignedVm, U, alingnednframes);
alignedVariance = nanmean(alignedVariance,1);

%% inner function
    function alignedVc = innerAlignFunc(Vc)
    alignedVc = NaN(size(Vc,1), segFrames(5), size(Vc,3)); %new Vc to capture max duration of each segment
    for iTrials = 1:length(cBhv.Rewarded)
        % get indices for current trial
        stimOn = cBhv.RawEvents.Trial{iTrials}.Events.Wire3High; %time of stimulus onset - measured from soundcard
        handleOn = [reshape(cBhv.RawEvents.Trial{iTrials}.States.WaitForAnimal1',1,[]) ...
            reshape(cBhv.RawEvents.Trial{iTrials}.States.WaitForAnimal2',1,[]) ...
            reshape(cBhv.RawEvents.Trial{iTrials}.States.WaitForAnimal3',1,[])];



        clear cIdx
        cIdx(1) = handleOn(find(handleOn == cBhv.RawEvents.Trial{iTrials}.States.WaitForCam(1))-1); %handle grab
        cIdx(2) = stimOn; %stim onset
        cIdx(3) = max(cat(2,cBhv.stimEvents{iTrials}{:})) + stimOn; %time of last stimulus event, start of delay
        cIdx(4) = cBhv.RawEvents.Trial{iTrials}.States.MoveSpout(1); %spouts in - response period
        cIdx = floor((cIdx - cIdx(1)) * SRATE)  + HANDLEGRABFRAME; %convert to frames relative to handle grab. This is the last frame of each segment.
        cIdx(5) = nframes; %the values in cIdx now specify when an event occurred on this trial

        alignedVc(:, 1 : segFrames(1), iTrials) = Vc(:, cIdx(1) - segFrames(1) + 1 : cIdx(1), iTrials); % baseline

        alignedVc(:, segFrames(1) + 1 : segFrames(1) + (diff(cIdx(1:2))), iTrials) = Vc(:, cIdx(1) + 1 : cIdx(2), iTrials); %handle period
        alignedVc(:, segFrames(2) + 1 : segFrames(2) + (diff(cIdx(2:3))), iTrials) = Vc(:, cIdx(2) + 1 : cIdx(3), iTrials); %stimulus period

        maxDiff = min([segFrames(3) + (diff(cIdx(3:4))) segFrames(5)]) - segFrames(3); %maximal possible delay duration

        if cIdx(4) > cIdx(5) %if response occurs after truncation of imaging data
            alignedVc(:, segFrames(3) + 1 : segFrames(3) + (cIdx(5) - cIdx(3)), iTrials) = Vc(:, cIdx(3) + 1 : cIdx(5), iTrials);
        else
            alignedVc(:, segFrames(3) + 1 : segFrames(3) + maxDiff, iTrials) = Vc(:, cIdx(3) + 1 : cIdx(3) + maxDiff, iTrials); %delay period
        end

        if segFrames(4) + (diff(cIdx(4:5))) > segFrames(5)
            alignedVc(:, segFrames(4) + 1 : segFrames(5), iTrials) = Vc(:, cIdx(4) + 1 : cIdx(4) + (segFrames(5) - segFrames(4)), iTrials); %response period
        else
            alignedVc(:, segFrames(4) + 1 : segFrames(4) + (diff(cIdx(4:5))), iTrials) = Vc(:, cIdx(4) + 1 : cIdx(5), iTrials); %response period
        end

    end
    end
end