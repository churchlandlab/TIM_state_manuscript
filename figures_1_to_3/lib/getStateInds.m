function [inds, attendinds,biasinds,Y, postprobs_sorted, bhv] = getStateInds(cPath,Animal,Rec,method,glmPath,dualCase, onlyLeftChoice)
Paradigm = 'SpatialDisc';
addpath('C:\Data\churchland\ridgeModel\rateDisc');

cPath = [cPath filesep Animal filesep Paradigm filesep Rec filesep]; %Widefield data path
bhvFile = dir([cPath filesep Animal '_' Paradigm '*.mat']);
load([cPath bhvFile(1).name],'SessionData'); %load behavior data
bhv = SessionData;clear SessionData;
load(glmPath,'posterior_probs','model_training_sessions','state_label_indices','mouse','masks'); %load behavior data
model_training_sessions = num2cell(model_training_sessions,2); %convert to a cell for ease
model_training_sessions = strtrim(model_training_sessions);
mouse = cellstr(mouse);

sessionind = find(strcmp(model_training_sessions,Rec) & strcmp(mouse, Animal));%find the index of the session we want to pull latent states for

postprob_nonan = posterior_probs{sessionind}; %grab the proper session
postprobs_withnan = NaN(length(state_label_indices) ,length(masks{sessionind}));
postprobs_withnan(:,masks{sessionind}) = postprob_nonan';

postprobs_sorted = postprobs_withnan(state_label_indices,:); %permute the states so theyre in the correct indices

if onlyLeftChoice
    useIdx = bhv.ResponseSide == 1; % 1 is left, 2 is right
else
    useIdx = ~isnan(bhv.ResponseSide); %only use performed trials
end

if strcmp(method,'max')
    [~,state1hot] = max(postprobs_sorted,[],1);
elseif strcmp(method,'cutoff')
    binary = postprobs_sorted' > .8; %get indices with P(state) > .8
    for i = 1:size(binary,1)
        temp = binary(i,:) == 1;
        if sum(temp) == 0
            state1hot(i) = 0; % if no states are P > .8, set to 0 to indicate no state
        else
            state1hot(i) = find(temp);
        end
    end
    Atemp = state1hot == 1;
    Btemp = state1hot == 2 | state1hot == 3;
    useIdx = useIdx & (Atemp | Btemp); %and only use trials where P of ANY state was > .8
end

if strcmp(dualCase,'choice')
    inds = find(rateDisc_equalizeTrials(useIdx, state1hot == 1, bhv.ResponseSide == 1, inf, true)); %equalize state AND L/R choices
elseif strcmp(dualCase,'reward')
    inds = find(rateDisc_equalizeTrials(useIdx, state1hot == 1, bhv.Rewarded == 1, inf, true));  %equalize to state AND rewarded vs unrewarded
elseif strcmp(dualCase,'none')
    inds = find(rateDisc_equalizeTrials(useIdx, state1hot == 1, [], inf, []));  %equalize to state only
elseif strcmp(dualCase,'choice_equal')
    inds = find(rateDisc_equalizeTrials(useIdx, state1hot == 1, bhv.ResponseSide == 1, inf, false)); %equalize state AND L/R choices
else
    error('Need to input a valid dualCase parameter')
end

attendinds = inds(state1hot(inds) == 1);
biasinds = inds(state1hot(inds) ~= 1);
Y = state1hot(inds) == 1;
end