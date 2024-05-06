function [DLCwavelet] = MovementFeatures(Lateral_allFrames, Bottom_allFrames, aligned_FrameTime, raw_data, HMM_states, outlier_filter)

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

figure('Name', 'Epoch time length distribution');
subplot(2,2,1);
hold on
histogram(length_stim_L);
histogram(length_stim_B);
title('Stimulus');
xlabel('Frame Number');
ylabel('Trial Number');
legend({'Lateral', 'Bottom'});
disp(['The biggest mis-match between lateral and bottom cameras for the stimulus window is ', num2str(max(abs(length_stim_L - length_stim_B)))]);



length_delay_L = aligned_FrameTime.spoutsIn(1,:) - aligned_FrameTime.stimOff(1,:);
length_delay_B = aligned_FrameTime.spoutsIn(2,:) - aligned_FrameTime.stimOff(2,:);

subplot(2,2,2);
hold on
histogram(length_delay_L);
histogram(length_delay_B);
title('Delay');
xlabel('Frame Number');
ylabel('Trial Number');
disp(['The biggest mis-match between lateral and bottom cameras for the delay window is ', num2str(max(abs(length_delay_L - length_delay_B)))]);



length_response_L = cellfun(@(x)(size(x,1)),Lateral_allFrames.frames,'UniformOutput',true) - aligned_FrameTime.spoutsIn(1,:);
length_response_B = cellfun(@(x)(size(x,1)),Bottom_allFrames.frames,'UniformOutput',true) - aligned_FrameTime.spoutsIn(2,:);

subplot(2,2,3);
hold on
histogram(length_response_L);
histogram(length_response_B);
title('Response');
xlabel('Frame Number');
ylabel('Trial Number');
disp(['The biggest mis-match between lateral and bottom cameras for the response window is ', num2str(max(abs(length_response_L - length_response_B)))]);



length_all_L = 15 + length_stim_L + length_delay_L + length_response_L;
length_all_B = 15 + length_stim_B + length_delay_B + length_response_B;

a = [length_all_L,length_all_B];
subplot(2,2,4);
hold on
histogram(length_all_L, [min(a):2:max(a)]);
histogram(length_all_B, [min(a):2:max(a)]);
title('All Epochs');
xlabel('Frame Number');
ylabel('Trial Number');
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


figure('Name', 'Averaged motion energy');
hold on
a = stdshade(motionEnergy_engaged(:,1:110), 0.4, [0 0.4470 0.7410]);
b = stdshade(motionEnergy_disengaged(:,1:110), 0.4, [0.8500 0.3250 0.0980]);
legend([a, b], {'Engaged', 'Disengaged'});
xlabel('Frame from Baseline');
set(gca,'box','off');
set(gca,'tickdir','out');
hold off



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



[~, p_value] = ttest2(motionA, motionB);
fig_energy = figure('Name', ['Motion energy distribution across states (stimulus+delay), ', mousename]);
subplot(1,2,1);
a = histogram(motionA, 0:0.02:1.5, 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceAlpha', 0.7);
hold on
b = histogram(motionB, 0:0.02:1.5, 'FaceColor', [0 0.4470 0.7410], 'EdgeColor', 'none', 'FaceAlpha', 0.8);
title(['All DLC labels, ', 'p value = ', num2str(p_value)]);

yylim = ylim .* 1.1;
line([mean(motionA) mean(motionA)], yylim, 'Color', 'red');
line([mean(motionB) mean(motionB)], yylim, 'Color', [0 0.4470 0.7410]);
ylim(yylim);

axis square
xlabel('Normalized Motion Energy');
ylabel('Number of Trials');

legend([a, b], {'Engaged', 'Disengaged'});

set(gca,'box','off');
set(gca,'tickdir','out');

hold off
clear a b




%% 4th, look into the frequency domain to test the temporarily consistent movement features (using the first 10 PCs of DLC positions)
% Here's something: doing one PCA on both engaged and disengaged trials will make the results comparable across states. 
% However, doing PCA separately may better capture the differences in frequency domain.
% Only use the stimulus + delay windows
spoutIn_all = mean(aligned_FrameTime.spoutsIn, 1) - mean(aligned_FrameTime.stimOn, 1) + 1;
longest_trial = max(floor(spoutIn_all));

DLC = areaNormalization(DLC);
DLC_forPCA = nan(size(DLC, 1), longest_trial, size(DLC, 3));

for i = 1 : num_trial
    t = floor(mean(aligned_FrameTime.spoutsIn(:, i), 1) - mean(aligned_FrameTime.stimOn(:, i), 1));
    DLC_forPCA(i, 1:t, :) = DLC(i, 16:t+15, :);
    
end


a = size(DLC_forPCA);
ttt = reshape(DLC_forPCA, [], a(3));

[~,PC,~,~,explained,~] = pca(ttt, 'Rows', 'complete');
PC = reshape(PC, a);


% a quick test on if the PCs carry important task information
PC_L = squeeze(PC(idx_L, 30, :));
PC_R = squeeze(PC(idx_R, 30, :));

figure('Name', ['PCA and choice, ', mousename]);
hold on
a = scatter3(PC_L(:,1), PC_L(:,2), PC_L(:,3), 10, 'filled');
b = scatter3(PC_R(:,1), PC_R(:,2), PC_R(:,3), 10, 'filled');
xlabel('PC 1');
ylabel('PC 2');
zlabel('PC 3');
legend([a, b], {'Left Choice', 'Right Choice'});
rotate3d on

hold off
clear a b




PC_engaged = PC(idx_engaged, :, 1:9);
PC_disengaged = PC(idx_disengaged, :, 1:9);




a = abs(diff(PC_engaged, 1, 2));
PCenergyA = nanmean(reshape(a, size(a, 1), []), 2);
b = abs(diff(PC_disengaged, 1, 2));
PCenergyB = nanmean(reshape(b, size(b, 1), []), 2);

if length(PCenergyA) > length(PCenergyB)  
    e = randsample(length(PCenergyA),length(PCenergyB));
    PCenergyA = PCenergyA(e);  
else
    e = randsample(length(PCenergyB), length(PCenergyA));
    PCenergyB = PCenergyB(e);
end
clear a b e


figure(fig_energy);
subplot(1,2,2);
hold on
histogram(PCenergyA, 0:0.05:5, 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceAlpha', 0.7);
hold on
histogram(PCenergyB, 0:0.05:5, 'FaceColor', [0 0.4470 0.7410], 'EdgeColor', 'none', 'FaceAlpha', 0.8);
title(['Fisrt 10 PCs, ', 'p value = ', num2str(p_value)]);

yylim = ylim .* 1.1;
line([mean(motionA) mean(motionA)], yylim, 'Color', 'red');
line([mean(motionB) mean(motionB)], yylim, 'Color', [0 0.4470 0.7410]);
ylim(yylim);

axis square
xlabel('Sum of Diff across Frames');
ylabel('Number of Trials');

set(gca,'box','off');
set(gca,'tickdir','out');

hold off
clear a b





figure('Name', ['PCA and states, ', mousename]);
hold on
a = scatter3(PC_engaged(:,30,1), PC_engaged(:,30,2), PC_engaged(:,30,3), 10, 'filled');
b = scatter3(PC_disengaged(:,30,1), PC_disengaged(:,30,2), PC_disengaged(:,30,3), 10, 'filled');
xlabel('PC 1');
ylabel('PC 2');
zlabel('PC 3');
legend([a, b], {'Engaged', 'Disengaged'});
rotate3d on

hold off
clear a b




figure('Name', ['DLC PCA results, ', mousename]);
for i = 1 : size(PC_engaged, 3)
    subplot(3,3,i);
    hold on
    title(['PC ', num2str(i)]);
    a = stdshade(PC_engaged(:,:, i), 0.4, [0 0.4470 0.7410]);
    b = stdshade(PC_disengaged(:,:, i), 0.4, [0.8500 0.3250 0.0980]);
    legend([a, b], {'Engaged', 'Disengaged'});
    hold off
end




sampleRate = 30;    % (Hz), Sampling rate of widefield imaging
windowSize = 1;     % (s), the width of convolution gaussian 
FWHM = 0.3;         % (s), Full width at half maximum, determining the width of gaussian 
min_freq = 0;       % (Hz), minimal frequency for time-frequency plot
max_freq = 5;      % (Hz), maximal frequency for time-frequency plot
num_freq = 10;      % number of frequency in count



engagedPCs_across_frequency = [];
colorRange = nan(9, 2); % This variable is to make sure the plots of two states are at the same scale.

figure('Name', ['DLC wavelet, engaged trials, ', mousename]);
for i = 1 : size(PC_engaged, 3)
    subplot(3,3,i);
    title(['PC ', num2str(i)]);
    % Please notice, the PC processed by wavelet decomposition is already smoothed by the Gaussian used in convolution
    [temporalFeatures_enaged, colorRange(i, :)] = Wavelet(PC_engaged(:,:,i)', sampleRate, windowSize, FWHM, min_freq, max_freq, num_freq, []);
    
    if i == 1
        engagedPCs_across_frequency = temporalFeatures_enaged;
    else
        engagedPCs_across_frequency = cat(1, engagedPCs_across_frequency, temporalFeatures_enaged);
    end

end



disengagedPCs_across_frequency = [];

figure('Name', ['DLC wavelet, disengaged trials, ', mousename]);
for i = 1 : size(PC_engaged, 3)
    subplot(3,3,i);
    title(['PC ', num2str(i)]);
    [temporalFeatures_disenaged, ~] = Wavelet(PC_disengaged(:,:,i)', sampleRate, windowSize, FWHM, min_freq, max_freq, num_freq, colorRange(i,:));
    
    if i == 1
        disengagedPCs_across_frequency = temporalFeatures_disenaged;
    else
        disengagedPCs_across_frequency = cat(1, disengagedPCs_across_frequency, temporalFeatures_disenaged);
    end
end



x = squeeze(nanmean(engagedPCs_across_frequency, 2));
y = squeeze(nanmean(disengagedPCs_across_frequency, 2));
DLCwavelet.PowerAcrossFrequency_engaged = x;
DLCwavelet.PowerAcrossFrequency_disengaged = y;


figure('Name', ['DLC wavelet power across frequencies, ', mousename]);
for i = 1 : size(PC_engaged, 3)
    subplot(3,3,i);
    hold on
    range = ((i-1)*num_freq:(i-1)*num_freq+(num_freq - 1)) + 1;
    
    a = stdshade(x(range,:)', 0.4, [0 0.4470 0.7410]);
    b = stdshade(y(range,:)', 0.4, [0.8500 0.3250 0.0980]);
    
    xlim([0 num_freq]);
    xticks([1,num_freq]);
    xticklabels({num2str(min_freq), num2str(max_freq)});
    
    xlabel('Frequency');
    ylabel('Power');
    if i == 1
        legend([a, b], {'Engaged', 'Disengaged'});
    end
    
    set(gca,'box','off');
    set(gca,'tickdir','out');
    hold off
end



x = engagedPCs_across_frequency;
y = disengagedPCs_across_frequency;
figure('Name', ['DLC wavelet power across time, ', mousename]);
for i = 1 : size(PC_engaged, 3)
    subplot(3,3,i);
    hold on
    range = ((i-1)*num_freq:(i-1)*num_freq+(num_freq - 1)) + 1;
    
    xx = squeeze(nanmean(x(range,:, :), 1));
    yy = squeeze(nanmean(y(range,:, :), 1));
    
    a = stdshade(xx', 0.4, [0 0.4470 0.7410]);
    b = stdshade(yy', 0.4, [0.8500 0.3250 0.0980]);
    
    
    xlabel('Frames');
    ylabel('Power');
    if i == 1
        legend([a, b], {'Engaged', 'Disengaged'});
    end
    
    set(gca,'box','off');
    set(gca,'tickdir','out');
    hold off
end




x = engagedPCs_across_frequency;
y = disengagedPCs_across_frequency;

x = squeeze(nanmean(nanmean(x, 2), 1));
y = squeeze(nanmean(nanmean(y, 2), 1));
[~, p_value] = ttest2(x, y);

figure('Name', ['DLC wavelet power across states, ', mousename]);
a = histogram(x, 0:0.03:max([x(:); y(:)]), 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceAlpha', 0.7);
hold on
b = histogram(y, 0:0.03:max([x(:); y(:)]), 'FaceColor', [0 0.4470 0.7410], 'EdgeColor', 'none', 'FaceAlpha', 0.8);
title(['p value = ', num2str(p_value)]);

yylim = ylim .* 1.1;
line([mean(x) mean(x)], yylim, 'Color', 'red');
line([mean(y) mean(y)], yylim, 'Color', [0 0.4470 0.7410]);
ylim(yylim);

axis square
xlabel('Averaged Wavelet Power');
ylabel('Number of Trials');

legend([a, b], {'Engaged', 'Disengaged'});

set(gca,'box','off');
set(gca,'tickdir','out');

hold off
clear a b






%% 5th, t-SNE on the wavelet decomposition results
engagedPCs_across_frequency = engagedPCs_across_frequency(:, 1:10:end, :);
disengagedPCs_across_frequency = disengagedPCs_across_frequency(:, 1:10:end, :);

allPCs_across_frequency = cat(3, engagedPCs_across_frequency, disengagedPCs_across_frequency);
index = size(engagedPCs_across_frequency, 3);

allPCs_across_frequency = permute(allPCs_across_frequency, [3 1 2]);
allPCs_across_frequency = reshape(allPCs_across_frequency, size(allPCs_across_frequency,1), []);



if size(allPCs_across_frequency, 2) > size(allPCs_across_frequency, 1)
    msgfig = msgbox('More dimensions than samples in t-SNE!','Error','modal');
    uiwait(msgfig);
    return
end


Y = tsne(allPCs_across_frequency);

engaged_tSNE = Y(1:index,:);
disengaged_tSNE = Y(index+1:end,:);
DLCwavelet.engaged_tSNE =  engaged_tSNE;
DLCwavelet.disengaged_tSNE = disengaged_tSNE;


% the next part is to calculate the spatial density of t-SNE results  
bin_x = (max(Y(:,1)) - min(Y(:,1))) / 30;
bin_y = (max(Y(:,2)) - min(Y(:,2))) / 30;


[N_engaged, ~] = hist3(engaged_tSNE, 'Ctrs',{min(Y(:,1)):bin_x:max(Y(:,1)) min(Y(:,2)):bin_y:max(Y(:,2))});
[N_disengaged, ~] = hist3(disengaged_tSNE, 'Ctrs',{min(Y(:,1)):bin_x:max(Y(:,1)) min(Y(:,2)):bin_y:max(Y(:,2))});

N_engaged = (N_engaged ./ index) .* 100;
N_disengaged = (N_disengaged ./ size(disengagedPCs_across_frequency, 3)) .* 100;




figure('Name', ['t-SNE result of DLC trajectories, ', mousename]); 
subplot(1,4,1); 
scatter(engaged_tSNE(:, 1), engaged_tSNE(:, 2), 10, 'red', 'filled');
hold on
scatter(disengaged_tSNE(:, 1), disengaged_tSNE(:, 2), 10, [0 0.4470 0.7410], 'filled');
xlabel('t-SNE axis 1');
ylabel('t-SNE axis 2');
legend({'Engaged', 'Disengaged'});
set(gca,'box','off');
set(gca,'tickdir','out');
hold off


a = imgaussfilt(flipud(N_engaged'), 1.5);
b = imgaussfilt(flipud(N_disengaged'), 1.5);
minimal = min([a(:); b(:)]);
maximal = max([a(:); b(:)]);

subplot(1,4,2);
imagesc(a, [minimal maximal]);
title(['Engaged, entropy=', num2str(entropy(a))]);
set(gca,'box','off');
set(gca,'tickdir','out');
hold off


subplot(1,4,3);  
imagesc(b, [minimal maximal]);
title(['Disngaged, entropy=', num2str(entropy(b))]);
set(gca,'box','off');
set(gca,'tickdir','out');
hold off

clear minimal maximal



% It seems that the t-SNE result of wavelet decomposition isn't very stable, we need to make some shuffle controls.
realDiff = mean(reshape(abs(a - b), 1, []));
shuffleDiff = nan(1, 1000);


for i = 1 : 1000
    
    allPCs_across_frequency = allPCs_across_frequency(randperm(size(allPCs_across_frequency, 1)), :);
    Y = tsne(allPCs_across_frequency);
    
    
    bin_x = (max(Y(:,1)) - min(Y(:,1))) / 30;
    bin_y = (max(Y(:,2)) - min(Y(:,2))) / 30;
    
    [N_engaged, ~] = hist3(Y(1:index,:), 'Ctrs',{min(Y(:,1)):bin_x:max(Y(:,1)) min(Y(:,2)):bin_y:max(Y(:,2))});
    [N_disengaged, ~] = hist3(Y(index+1:end,:), 'Ctrs',{min(Y(:,1)):bin_x:max(Y(:,1)) min(Y(:,2)):bin_y:max(Y(:,2))});
    
    N_engaged = (N_engaged ./ index) .* 100;
    N_disengaged = (N_disengaged ./ size(disengagedPCs_across_frequency, 3)) .* 100;
    
    
    a = imgaussfilt(flipud(N_engaged'), 1.5);
    b = imgaussfilt(flipud(N_disengaged'), 1.5);
    
    shuffleDiff(i) = mean(reshape(abs(a - b), 1, []));
    clear a b bin_x bin_y Y N_engaged N_disengaged
    
end




subplot(1,4,4);
hold on
h = histogram(shuffleDiff, 'FaceColor', [0 0.4470 0.7410], 'EdgeColor', 'none', 'FaceAlpha', 0.8);
aaa = line([realDiff realDiff], ylim, 'Color','red');

if kstest(normalize(shuffleDiff)) == 1 % One-sample Kolmogorov-Smirnov test
    title('The shuffle distribution is NOT Gaussian');
else
    p = 1 - normcdf(realDiff, mean(shuffleDiff), std(shuffleDiff));
    title(['The control distribution is Gaussian, p value = ', num2str(p)]);
end




x = min(shuffleDiff) : 0.0001 : max(shuffleDiff);
y = normpdf(x,mean(shuffleDiff),std(shuffleDiff));

ttt = max(h.Values) / max(y);

plot(x, y.*ttt,'Color', [0.6 0.6 0.6], 'LineWidth', 2);

legend(aaa, {'True Diff across States'});
xlabel('t-SNE Density Map Diff');
ylabel('Number of Shuffles');
set(gca,'box','off');
set(gca,'tickdir','out');
hold off

save(['X:\Chaoqun\Papers&Talks\TIV-state Manuscript\Wavelet_fig\', mouse_name, '_waveletResults.mat'], 'engaged_tSNE', 'disengaged_tSNE', 'realDiff', 'shuffleDiff', "-mat");


%% Result saving
DLCwavelet.shuffleControl = shuffleDiff;
DLCwavelet.realDiff = realDiff;



end






function [temporalFeatures, colorRange] = Wavelet(data, sampleRate, windowSize, FWHM, min_freq, max_freq, num_freq, colorRange)


if max_freq > (sampleRate/2)
    msgfig = msgbox('The max frequency must NOT be higher than 1/2 of sampling rate!','Error','modal');
    uiwait(msgfig);
    return
end


    
% 1st, setting up the time length of convolution window
time = (0:windowSize*sampleRate-1)/sampleRate;
time = time - mean(time); % best practice is to have time=0 at the center of the wavelet


dataR = reshape(data,1,[]);
dataR(isnan(dataR)) = nanmean(dataR);   % Please notice, all NaNs in the PCA results get replaced with mean values in this step, 
                                        % so the Wavelet output will have no NaNs.

%2nd, convolution
% Step 1: N's of convolution
ndata = length(dataR); 
nkern = length(time);
nConv = ndata + nkern;
halfK = floor(nkern/2);

% Step 2: FFTs
dataX = fft(dataR,nConv);


frex = linspace(min_freq,max_freq,num_freq);

% initialize TF matrix
timevec = [1 : size(data, 1)] ./ sampleRate;

tf = zeros(num_freq,length(timevec));


temporalFeatures = nan([num_freq, size(data)]);

for fi=1 : num_freq
    
    % create wavelet
    cmw  = exp(1i*2*pi*frex(fi)*time) .* ...
           exp( -4*log(2)*time.^2 / FWHM^2 );
    
    cmwX = fft(cmw,nConv);
    cmwX = cmwX./max(abs(cmwX));    % normalize the wavelet
    
    % FDplot(cmw, cmwX, dataX, time, sampleRate, nConv);   % this FDplot functions plots the gaussian, data, and filtered data in the frequency domain. 
                % You can use this function to check if the result makes sense.
    
    % the rest of convolution
    as = ifft( dataX.*cmwX );
    as = as(halfK+1:end-halfK);
%     as(end + 1) = mean(as);
    as = reshape(as,size(data));
    
    % extract power
    aspow = abs(as).^2;
    % aspow = abs(as);
    temporalFeatures(fi, :, :) = aspow;
    
    % average over trials and put in matrix
    tf(fi,:) = mean(aspow,2);
    
end




hold on
ttt = tf(:);

if isempty(colorRange)
    colorRange = [min(ttt), max(ttt)];
end

contourf(timevec,frex,tf,40,'linecolor','none')
set(gca,'ylim',[min_freq max_freq], 'clim', colorRange);
xlabel('Time (s)'), ylabel('Frequency (Hz)');
hold off


end

