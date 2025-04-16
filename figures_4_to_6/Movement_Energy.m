function Movement_Energy(Lateral_allFrames, Bottom_allFrames, aligned_FrameTime, raw_data, HMM_states, outlier_filter)

global mouse_name Laterl_labels Bottom_labels
mousename = mouse_name; 
DLC_labels = [Laterl_labels, Bottom_labels];

num_trial = size(aligned_FrameTime.stimOn, 2);
num_label_L = length(Laterl_labels);
num_label_B = length(Bottom_labels);


if isempty(HMM_states)
    xba = smoothdata(raw_data.Rewarded, 'gaussian', 50);
else
    xba = smoothdata(HMM_states, 'gaussian', 50);
end



%% 1st, showing the distribution of the timelength of different task epochs. (baseline is always 15 frames.)
length_stim_L = aligned_FrameTime.stimOff(1,:) - aligned_FrameTime.stimOn(1,:);
length_stim_B = aligned_FrameTime.stimOff(2,:) - aligned_FrameTime.stimOn(2,:);


length_delay_L = aligned_FrameTime.spoutsIn(1,:) - aligned_FrameTime.stimOff(1,:);
length_delay_B = aligned_FrameTime.spoutsIn(2,:) - aligned_FrameTime.stimOff(2,:);


length_response_L = cellfun(@(x)(size(x,1)),Lateral_allFrames.frames,'UniformOutput',true) - aligned_FrameTime.spoutsIn(1,:);
length_response_B = cellfun(@(x)(size(x,1)),Bottom_allFrames.frames,'UniformOutput',true) - aligned_FrameTime.spoutsIn(2,:);


length_all_L = 15 + length_stim_L + length_delay_L + length_response_L;
length_all_B = 15 + length_stim_B + length_delay_B + length_response_B;

disp(['The biggest mis-match between lateral and bottom cameras for the entire trial is ', num2str(max(abs(length_all_L - length_all_B)))]);




%% 2nd, align all trials to Baseline.
max_frameNum = max([length_all_L, length_all_B]);

DLC_L = nan(num_trial, max_frameNum, num_label_L*2);
DLC_B = nan(num_trial, max_frameNum, num_label_B*2);

for i = 1 : num_trial
    ttt = Lateral_allFrames.frames{1, i};
    ttt(:, 3:3:end) = [];
    DLC_L(i, 1:length_all_L(i), :) = ttt(aligned_FrameTime.stimOn(1, i)-14:end, :);
    
    ttt2 = Bottom_allFrames.frames{1, i};
    ttt2(:, 3:3:end) = [];
    DLC_B(i, 1:length_all_B(i), :) = ttt2(aligned_FrameTime.stimOn(2, i)-14:end, :);
end
clear i ttt ttt2


DLC = cat(3, DLC_L, DLC_B);





%% 3rd, motion energy sorted based on peak time before response window & grouped based on states
motionEnergy = diff(DLC, 1, 2);
motionEnergy = sqrt(motionEnergy(:, :, 1:2:end).^2 + motionEnergy(:, :, 2:2:end).^2);

motionEnergy = squeeze(nanmean(motionEnergy, 3));


% Applying the outliner filter here
if outlier_filter > 0
    all_motion = nanmean(motionEnergy, 2);
    mu = mean(all_motion);
    sigma = std(all_motion);
    
    idx_outliner = [find(all_motion < (mu - 3*sigma)); find(all_motion > (mu + 3*sigma))];
    
    
    ttt = motionEnergy;
    ttt(idx_outliner, :) = [];
    motionEnergy(idx_outliner, :) = repelem(nanmean(ttt, 1), length(idx_outliner), 1);
    
    ttt = DLC;
    ttt(idx_outliner, :, :) = [];
    DLC(idx_outliner, :, :) = repelem(nanmean(ttt, 1), length(idx_outliner), 1, 1);
end
clear ttt mu sigma all_motion

size_a = size(motionEnergy);
B = reshape(motionEnergy, 1, []);
B = normalize(B, 'zscore');
motionEnergy = reshape(B, size_a);
clear size_a B



idx_L = raw_data.ResponseSide == 1;
idx_R = raw_data.ResponseSide == 2;

choicehistory = [0, raw_data.ResponseSide];
choicehistory(end) = [];
idx_Lpre = find(choicehistory == 1);
idx_Rpre = find(choicehistory == 2);



[B,I] = sort(xba);

idx_disengaged = I(1 : floor(length(I) / 5));   % Use the trials with top 20% HMM value as engaged trials.
idx_engaged = I(end - floor(length(I) / 5) : end); % Use the trials with bottom 20% HMM value as disengaged trials.


motionEnergy_disengaged = motionEnergy(idx_disengaged, :);
motionEnergy_engaged = motionEnergy(idx_engaged, :);
clear B I





spoutIn_engaged = mean(aligned_FrameTime.spoutsIn(:, idx_engaged), 1) - mean(aligned_FrameTime.stimOn(:, idx_engaged), 1) + 15;
spoutIn_disengaged = mean(aligned_FrameTime.spoutsIn(:, idx_disengaged), 1) - mean(aligned_FrameTime.stimOn(:, idx_disengaged), 1) + 15;



for_colorBar = sort(motionEnergy(:));
for_colorBar(isnan(for_colorBar)) = [];
a = median(motionEnergy_engaged(:),'omitnan') + ((median(motionEnergy_engaged(:),'omitnan') - min(for_colorBar)))*2.2;    % the median value may be negative.
for_colorBar = [min(for_colorBar)*1.15 a];


figure('Name', ['Motion energy sorted based on spoutIn time, ', mouse_name]);
subplot(1,2,1);
hold on
[B, I] = sort(spoutIn_engaged, 'descend');
X = motionEnergy_engaged(I, :);
imagesc(X(:, 1:120), for_colorBar);
line([15 15], ylim, 'Color', 'white', 'LineWidth', 2);
scatter(B, 1:size(X, 1), 4, 'filled', 'red');
xlim([1 120]);
ylim([1 size(X, 1)]);
title('Engaged');
xlabel('Frame from Baseline');
ylabel('Trials');
set(gca,'box','off');
set(gca,'tickdir','out');
colorbar()
hold off


subplot(1,2,2);
hold on
[B, I] = sort(spoutIn_disengaged, 'descend');
X = motionEnergy_disengaged(I, :);
imagesc(X(:, 1:120), for_colorBar);
line([15 15], ylim, 'Color', 'white', 'LineWidth', 2);
scatter(B, 1:size(X, 1), 4, 'filled', 'red');
xlim([1 120]);
ylim([1 size(X, 1)]);
title('Disengaged');
xlabel('Frame from Baseline');
ylabel('Trials');
set(gca,'box','off');
set(gca,'tickdir','out');
colorbar
hold off




motionA = nan(size(spoutIn_engaged));
for i = 1 : length(spoutIn_engaged)
    motionA(i) = nanmean(motionEnergy_engaged(i, 16:floor(spoutIn_engaged(i))));
end

motionB = nan(size(spoutIn_disengaged));
for i = 1 : length(spoutIn_disengaged)
    motionB(i) = nanmean(motionEnergy_disengaged(i, 16:floor(spoutIn_disengaged(i))));
end



if length(motionA) > length(motionB)
    
    a = randsample(length(motionA),length(motionB));
    motionA = motionA(a);
    
else
    a = randsample(length(motionB), length(motionA));
    motionB = motionB(a);
end
clear a i


end





