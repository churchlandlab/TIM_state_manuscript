function Movement_DLC_analysis(fPath, mousename)


global mouse_name   % this name is widely used in following functions
mouse_name = mousename;



%% 1st, loading DeepLabCut (DLC) movement tracking and aligned camera frame times.
[Lateral_allFrames, Bottom_allFrames] = loadDLC(fPath, mousename);


aligned_FrameTime = dir([fPath mousename '_SpatialDisc_*' 'cameraTimes.mat']);
aligned_FrameTime = load([fPath aligned_FrameTime.name]);  % 1st row for lateral camera, 2nd row for bottom one 
aligned_FrameTime = aligned_FrameTime.aligned_FrameTime; 




%% 2nd, pre-processing DLC tracking, GLM-HMM state, and other data
Behave_Date = dir([fPath mousename '_SpatialDisc_*' 'frameTimes']);
Behave_Date_file = dir([fPath Behave_Date.name '\' mousename '_SpatialDisc_*' 'Session1.mat']);

Behave_Date = load([Behave_Date_file.folder '\' Behave_Date_file.name]);
Behave_Date = Behave_Date.SessionData;



FilteredPupil = load([Behave_Date_file.folder '\' 'FilteredPupil.mat']); % Loading the pupil diameter file, which is named "FilteredPupil"

if isfield(FilteredPupil, 'all_pupil') == 1
    FilteredPupil = FilteredPupil.all_pupil; 
elseif isfield(FilteredPupil, 'FilteredPupil') == 1
    FilteredPupil = FilteredPupil.FilteredPupil;
end


  
[Lateral_allFrames, Bottom_allFrames, aligned_FrameTime, Behave_Date, FilteredPupil, delete_idx] = FrameProcessing...
    (fPath, Lateral_allFrames, Bottom_allFrames, aligned_FrameTime, Behave_Date, FilteredPupil); % pre-processing the frames




HMM_state = load(['X:\StateProjectCentralRepo\DLC_results\', mousename, '\HMM_state_', mousename, '.mat']);
HMM_state = HMM_state.HMM_state;
HMM_states = HMMstateProcessing(HMM_state, delete_idx, Behave_Date);



%% 3rd, movement analyses
outlier_filter = 5; % Removing the trials with DLC labels that are 5 std away from the mean position

Movement_Energy(Lateral_allFrames, Bottom_allFrames, aligned_FrameTime, Behave_Date, HMM_states, outlier_filter);




%% 4th, task-independent-movement (TIM) calculation 
[Outcome, DLC, Fitted, DLCEnergy, TIM, State_results, borders] = TIM_analyses(Lateral_allFrames, Bottom_allFrames, aligned_FrameTime, Behave_Date, 50, outlier_filter, HMM_states);




%% 5th, pupil analyses
fillPupil = PupilProcessing(FilteredPupil, aligned_FrameTime, Lateral_allFrames); % Processing the pupil diameter data.
Pupil_plot = Pupil_performance_state(fillPupil, Behave_Date, HMM_states, 50);    % The last number is the size of the smooth window




%% 6th, save the results for further analyses (combining animals together)
save(['C:\Users\churchland\Desktop\State-TIM paper double check\', mouse_name, '_shared.mat'], 'Outcome', 'DLC', 'Fitted', 'DLCEnergy', 'TIM', 'State_results', 'borders', "-mat");

end