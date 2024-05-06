%% This code is adapted based on the original Mice_Behavior_Analysis code to test if doing the regression with separate 
%% session data is better. 2 major changes were made:
%%      1) The TIM regression was done separately for each session.
%%      2) Adding cross-validation to the regression.

function [] = Regression_Test(fPath, mousename)

%% 1st, loading DeepLabCut output files and aligned camera frame times.
% the later one can be got by running 'Behavior_alignVideoFrames.m' from Simon
[Lateral_allFrames, Bottom_allFrames] = getLabelEachFrame(fPath, mousename);

% back to upper folder
downpath = regexp(fPath,'\','split');
uppath = [];
for a = 1 : (size(downpath,2) - 2)
    uppath = [uppath cell2mat(downpath(a)) '\'];
end


aligned_FrameTime = dir([uppath mousename '_SpatialDisc_*' 'cameraTimes.mat']);
aligned_FrameTime = load([uppath aligned_FrameTime.name]);  % 1st row for lateral camera, 2nd row for bottom one 
clear downpath a

global mouse_name   % this name is widely used in following functions
mouse_name = mousename;



%% 2nd, extacting movement information from the output of DLC
raw_data = dir([uppath mousename '_SpatialDisc_*' 'frameTimes']);
raw_data_file = dir([uppath raw_data.name '\' mousename '_SpatialDisc_*' 'Session1.mat']);
if isempty(raw_data_file)
    raw_data_file = dir([uppath raw_data.name '\' mousename '_SpatialDisc_*' 'Session2.mat']);
    if isempty(raw_data_file)
        raw_data_file = dir([uppath raw_data.name '\' mousename '_SpatialDisc_*' 'Session3.mat']);
    end
end
raw_data = load([raw_data_file.folder '\' raw_data_file.name]);
raw_data = raw_data.SessionData;


if contains(fPath,'Jan32','IgnoreCase',true)
    aligned_FrameTime = aligned_FrameTime.aligned_FrameTime; 
elseif contains(fPath,'Jan33','IgnoreCase',true) 
    aligned_FrameTime = aligned_FrameTime.aligned_FrameTime;
end






% Loading the pupil diameter file, which is named "FilteredPupil" by Simon
pupil_file = dir([raw_data_file.folder '\' 'FilteredPupil.mat']);
if ~isempty(pupil_file)
    FilteredPupil = load([raw_data_file.folder '\' 'FilteredPupil.mat']);
    
    if isfield(FilteredPupil, 'all_pupil') == 1
        FilteredPupil = FilteredPupil.all_pupil;
        
    elseif isfield(FilteredPupil, 'FilteredPupil') == 1
        FilteredPupil = FilteredPupil.FilteredPupil;
    end
    
    
else
    FilteredPupil = [];
end
V1_activity = [];
M2_activity = [];


   
[Lateral_allFrames, Bottom_allFrames, aligned_FrameTime, raw_data, ~, ~, FilteredPupil, delete_idx] = FrameProcessing_regressiontest...
    (fPath, Lateral_allFrames, Bottom_allFrames, aligned_FrameTime, raw_data, V1_activity, M2_activity, FilteredPupil); % pre-processing the frames

clear raw_data_file downpath






outlier_filter = 5; % The number here denotes the standard deviation
% HMM_states = [];

HMM_state = load(['X:\StateProjectCentralRepo\DLC_results\', mousename, '\HMM_state_', mousename, '.mat']);
HMM_state = HMM_state.HMM_state;

HMM_states = HMMstateProcessing(HMM_state, delete_idx, raw_data);


DLCwavelet = MovementFeatures(Lateral_allFrames, Bottom_allFrames, aligned_FrameTime, raw_data, HMM_states, outlier_filter);

% [Outcome, DLC, Fitted, DLCEnergy, TIM, Results, State_results, borders] = LinearRegression_SD_regressiontest(Lateral_allFrames, Bottom_allFrames, aligned_FrameTime, raw_data, 50, outlier_filter, HMM_states);
[Outcome, DLC, Fitted, DLCEnergy, TIM, Results, State_results, borders] = LinearRegression_SD_final(Lateral_allFrames, Bottom_allFrames, aligned_FrameTime, raw_data, 50, outlier_filter, HMM_states);



fillPupil = PupilProcessing(FilteredPupil, aligned_FrameTime, Lateral_allFrames); % Processing the pupil diameter data.
Pupil_plot = Pupil_performance_state(fillPupil, raw_data, HMM_states, 50);    % The last number is the size of the smooth window


save(['X:\Chaoqun\Papers&Talks\TIV-state Manuscript\Fig5, major_TIM_figure\', mouse_name, '_final.mat'], 'Outcome', 'DLC', 'Fitted', 'DLCEnergy', 'TIM', 'Results', 'State_results', 'borders', "-mat");
% save(['X:\Chaoqun\Papers&Talks\TIV-state Manuscript\Fig6, Pupil_fig\', mouse_name, '_Pupil_plot.mat'], 'Pupil_plot', "-mat");

end