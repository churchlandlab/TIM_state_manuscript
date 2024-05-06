% This code is designed for pre-processing mouse video frames before using them

function [Lateral_allFrames, Bottom_allFrames, aligned_FrameTime, raw_data, V1_activity, M1_activity, FilteredPupil, delete_idx] = FrameProcessing_regressiontest...
    (fPath, Lateral_allFrames, Bottom_allFrames, aligned_FrameTime, raw_data, V1_activity, M1_activity, FilteredPupil)

global mouse_name  
mousename = mouse_name;

%% 1st, Checking tiral numbers
a = length(Lateral_allFrames.frames);
b = length(Bottom_allFrames.frames);
c = length(aligned_FrameTime.cameraTimes);
if a~=b || a~=c || b~=c
    message = 'Please make sure the lateral camera file, bottom camera file, and frame alignment file have SAME trial numbers!';
    msgfig = msgbox(message,'Error','modal');
    uiwait(msgfig);
    return;
end
clear a b c


%% 2nd, Removing the two spouts labels in bottom videos
for a = 1 : size(Bottom_allFrames.frames,2)
    Bottom_allFrames.frames{1, a}(:, 31:36) = [];
    Bottom_allFrames.frames{1, a}(:, 22:24) = [];
    % 22~24 denote the locations and likelihood of the tongue, which is supposed to be highly correlated the stimulus/choice. 
    % Thus, here is a concern: The negative correlation between task-independent variance and performance 
    % merely comes from the movement of the tongue. Of course, the correlation between condition-independent variance 
    % and the movement of tongue is negative. 
    % (Because the tongue is invisiable in most of the frames, it may also be a major noise recourse)
    % But, actually, it doesn't matter much. Tested on Fez7, CTP7, mSM63.
end


%% 3rd, Removing trials missing event markers
if isfield(raw_data,'optoType') && sum(~isnan(raw_data.optoType)) ~= 0 % This part is only for optogenetic trials
    raw_data.optoArea = raw_data.optoType;
    
    if length(raw_data.Notes) - length(raw_data.optoType) == 1
        raw_data.Notes(end) = [];   % There is one more empty note in the un-merged rawdata file
    elseif length(raw_data.Notes) - length(raw_data.optoType) == 0
        % Doing nothing in this case
    else
        fprintf(2, 'The raw_data.Notes of this file has some problems! \n');
        return;
    end
    
    for a = 1 : length(raw_data.Notes)
        if strcmpi(raw_data.Notes{a}, 'frontal')
            area_label = -1; % -1 denotes doing opto stimulation in frontal area
            control_label = -100; % -100 denotes this is a control trial in frontal stimulation sessions
        elseif strcmpi(raw_data.Notes{a}, 'parietal')
            area_label = 1;  % 1 denotes doing opto stimulation in parietal area
            control_label = 100; % 100 denotes this is a control trial in parietal stimulation sessions
        end
        
        
        if ~isnan(raw_data.optoArea(a))
            raw_data.optoArea(a) = area_label;
        elseif isnan(raw_data.optoArea(a))
            raw_data.optoArea(a) = control_label;
        end
       
    end
    clear area_label control_label
    
end

  
original_OrderIdx.Animal = mousename;
if isfield(raw_data, 'SessionNames')
    original_OrderIdx.SessionNames = raw_data.SessionNames;
end
original_OrderIdx.OrderIdx = [];
for i = 1 : length(raw_data.nTrials)
    original_OrderIdx.OrderIdx = [original_OrderIdx.OrderIdx, 1 : (raw_data.nTrials(i))];
end
if length(original_OrderIdx.OrderIdx) ~= sum(raw_data.nTrials)
    disp('The indice of trials have some problems!');
end
% This variable "original_OrderIdx" is designed for decoding task-related signals from widefield imaging data.
% In terms of the huge size of imaging matrix, we only be able to select the trials we use 
% after knowing which trial has high TIV or low TIV.
% In this situation, becuse we deleted some trials during the analyses, 
% we need the orignial order indice of the trials in their own sessions.
% "original_OrderIdx" is used for recording the original indice, the session names, and the number of trials of each session.


idx = find(isnan(aligned_FrameTime.stimOn(1,:)));
idx = [idx, find(isnan(aligned_FrameTime.stimOff(1,:)))];
idx = [idx, find(isnan(aligned_FrameTime.handlesIn(1,:)))];
idx = [idx, find(isnan(aligned_FrameTime.spoutsIn(1,:)))];
idx = unique(idx);
j = cumsum(raw_data.nTrials);
j = [0,j];
for a = 1 : (length(j)-1)
    raw_data.nTrials(a) = raw_data.nTrials(a) - length(idx(j(a)<idx & idx<=j(a+1)));
end
 
Lateral_allFrames.frames(idx) = [];
Lateral_allFrames.filename(idx) = [];
Bottom_allFrames.frames(idx) = [];
Bottom_allFrames.filename(idx) = [];
if ~isempty(V1_activity) || ~isempty(M1_activity)
    V1_activity(:,:,:,idx) = [];
    M1_activity(:,:,:,idx) = [];
end

fn = fieldnames(aligned_FrameTime);
for k=1:numel(fn)
    aligned_FrameTime.(fn{k})(:,idx) = []; 
end

fn = fieldnames(raw_data);
for k=1:numel(fn)
    if length(raw_data.(fn{k})) == (length(Lateral_allFrames.frames) + length(idx))
        raw_data.(fn{k})(:,idx) = [];
    end
end

original_OrderIdx.OrderIdx(idx) = [];

if ~isempty(FilteredPupil)
    FilteredPupil.orgPupil(idx) = [];
    FilteredPupil.fillPupil(idx) = [];
%     FilteredPupil.sPupil(idx) = [];
%     FilteredPupil.fPupil(idx) = [];
end

delete_idx = idx;
clear i fn k j a idx





global Original_OrderIdx   % this variable is widely used in following analyses
Original_OrderIdx = original_OrderIdx;


%% 4th, Normalization
global Lateral_allFrames_raw Bottom_allFrames_raw       
% Keep the raw DLC results without normalization for calculating the motion energy & average movetion
Lateral_allFrames_raw = Lateral_allFrames;
Bottom_allFrames_raw = Bottom_allFrames;

label_num_lateral = size(Lateral_allFrames.frames{1,1},2)/3;
label_num_bottom = size(Bottom_allFrames.frames{1,1},2)/3;
for a = 1 : label_num_lateral
    normal_array_x = [];
    normal_array_y = [];
    trial_length = [];
    idx_session = cumsum(raw_data.nTrials);     % Because the camera may be located slightly differnet in different days, we need to remove the effect of it. 
    idx = 1;
    for b = 1 : size(Lateral_allFrames.frames,2)
        normal_array_x = [normal_array_x; Lateral_allFrames.frames{1,b}(:, (a*3-2))];
        normal_array_y = [normal_array_y; Lateral_allFrames.frames{1,b}(:, (a*3-1))];
        trial_length = [trial_length, size(Lateral_allFrames.frames{1,b},1)];
        
        if ~isempty(find(idx_session == b, 1))
            normal_array_x(idx:end) = normal_array_x(idx:end) - mean(normal_array_x(idx:end));
            normal_array_y(idx:end) = normal_array_y(idx:end) - mean(normal_array_y(idx:end));
            idx = length(normal_array_x) + 1;
        end
    end

    
    t = cumsum(trial_length);
    t = [0, t];
    for b = 1 : size(Lateral_allFrames.frames,2)
        Lateral_allFrames.frames{1,b}(:, (a*3-2)) = normal_array_x( t(b)+1 : t(b+1) );
        Lateral_allFrames.frames{1,b}(:, (a*3-1)) = normal_array_y( t(b)+1 : t(b+1) );
    end
end


for a = 1 : label_num_bottom
    normal_array_x = [];
    normal_array_y = [];
    trial_length = [];
    idx_session = cumsum(raw_data.nTrials);
    idx = 1;
    for b = 1 : size(Bottom_allFrames.frames,2)
        normal_array_x = [normal_array_x; Bottom_allFrames.frames{1,b}(:, (a*3-2))];
        normal_array_y = [normal_array_y; Bottom_allFrames.frames{1,b}(:, (a*3-1))];
        trial_length = [trial_length, size(Bottom_allFrames.frames{1,b},1)];
        
        if ~isempty(find(idx_session == b, 1))
            normal_array_x(idx:end) = normal_array_x(idx:end) - mean(normal_array_x(idx:end));
            normal_array_y(idx:end) = normal_array_y(idx:end) - mean(normal_array_y(idx:end));
            idx = length(normal_array_x) + 1;
        end
    end

    
    t = cumsum(trial_length);
    t = [0, t];
    for b = 1 : size(Bottom_allFrames.frames,2)
        Bottom_allFrames.frames{1,b}(:, (a*3-2)) = normal_array_x( t(b)+1 : t(b+1) );
        Bottom_allFrames.frames{1,b}(:, (a*3-1)) = normal_array_y( t(b)+1 : t(b+1) );
    end
end


%% 5th, Getting the names of body parts
lateral_path = [fPath 'Lateral\'];
bottom_path = [fPath 'Bottom\'];

lateral_videos = dir([lateral_path mousename '_SpatialDisc_*.csv']);
bottom_videos = dir([bottom_path mousename '_SpatialDisc_*.csv']);


lateral_labelNames = Import_LabelNames([lateral_path lateral_videos(20).name], 2, 3);
lateral_labelNames(:, 43:end) = [];
bottom_labelNames = Import_LabelNames([bottom_path bottom_videos(20).name], 2, 3);  % just use the 20th video's label names
bottom_labelNames(:, 31:36) = [];
bottom_labelNames(:, 22:24) = [];
clear temp lateral_path bottom_path lateral_videos bottom_videos

global Laterl_labels Bottom_labels
Laterl_labels = unique(lateral_labelNames(1,:), 'stable');
Bottom_labels = unique(bottom_labelNames(1,:), 'stable');


end