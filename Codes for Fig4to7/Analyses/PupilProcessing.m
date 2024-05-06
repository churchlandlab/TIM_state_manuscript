%% This code is designed to find the different event windows of each trial in the FilteredPupil file according to aligned_FrameTime

function [fillPupil] = PupilProcessing(FilteredPupil, aligned_FrameTime, Lateral_allFrames)

fillPupil = {};
if size(aligned_FrameTime.stimOn, 2) == length(FilteredPupil.fillPupil)
    
    for i = 1 : length(FilteredPupil.fillPupil)
        
        if length(FilteredPupil.fillPupil{1, i}) == 500 % This is the max frame number set by Simon's Behavior_computePupilVars code
                                                        % if the orignial video has more frames, the code will only keep the last 500
                                                        % frames to get the pupil diameter.
            original_length = size(Lateral_allFrames.frames{1, i}, 1);
            diff = original_length - 500;
            aligned_FrameTime.stimOn(1, i) = aligned_FrameTime.stimOn(1, i) - diff;
            aligned_FrameTime.spoutsIn(1, i) = aligned_FrameTime.spoutsIn(1, i) - diff;
            Lateral_allFrames.frames{1, i}(1 : diff, :) = [];                                               
        end
        
        if abs(length(FilteredPupil.fillPupil{1, i}) - size(Lateral_allFrames.frames{1, i}, 1)) < 3 % If the difference between the alignments is less than 3 frames, it should be fine.
            Idx_stimOn = aligned_FrameTime.stimOn(1, i);
            Idx_spoutsIn = aligned_FrameTime.spoutsIn(1, i);
            Idx_delay = aligned_FrameTime.spoutsIn(1, i) - 15;
            Idx_handlesIn = aligned_FrameTime.handlesIn(1, i);
            
            fillPupil.Baseline{1, i} = FilteredPupil.fillPupil{1,i}(Idx_handlesIn - 15 : Idx_handlesIn - 1);
            
            % Please notice the time lengths of stimulus and delay windows are variable in many sessions
            fillPupil.StimuWindow{1, i} = FilteredPupil.fillPupil{1,i}(Idx_stimOn : Idx_stimOn + 29);
            fillPupil.DelayWindow{1, i} = FilteredPupil.fillPupil{1,i}(Idx_delay : Idx_delay + 14);
            
            fillPupil.RewardWindow{1, i} = FilteredPupil.fillPupil{1,i}(Idx_spoutsIn : end);
            
        else
            message = ['FilteredPupil and Lateral_allFrames have different frame numbers at trial ', num2str(i), '.'];
            msgfig = msgbox(message,'Error','modal');
            uiwait(msgfig);
            return
        end
    end
    
else
    msgfig = msgbox('FilteredPupil and aligned_FrameTime have different trial numbers!','Error','modal');
    uiwait(msgfig);
    return
end


end