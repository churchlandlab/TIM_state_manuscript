function [aligned_labels] = alignDLCToBehavior(labels, alignmentFrames, preFrames, postFrames)
aligned_labels = NaN(size(labels,2), size(labels{1},2), preFrames + postFrames);


for i = 1:length(labels)
    trialdata = labels{i};
    eventframe = alignmentFrames(i);
    if isnan(eventframe)
        continue
    end
    
    frames2grab =  eventframe - preFrames : eventframe + postFrames - 1;
    if sum(frames2grab < 1) > 0
        fprintf('skipping DLC alignment for this trial\n');
        continue
    end
    aligned_labels(i,:,:) = trialdata(frames2grab,:)';
end
end