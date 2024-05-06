mSM63 = load('X:\Chaoqun\Papers&Talks\TIV-state Manuscript\Fig6, Pupil_fig\mSM63_Pupil_plot.mat');
mSM64 = load('X:\Chaoqun\Papers&Talks\TIV-state Manuscript\Fig6, Pupil_fig\mSM64_Pupil_plot.mat');
mSM65 = load('X:\Chaoqun\Papers&Talks\TIV-state Manuscript\Fig6, Pupil_fig\mSM65_Pupil_plot.mat');
mSM66 = load('X:\Chaoqun\Papers&Talks\TIV-state Manuscript\Fig6, Pupil_fig\mSM66_Pupil_plot.mat');

mSM63 = mSM63.Pupil_plot;
mSM64 = mSM64.Pupil_plot;
mSM65 = mSM65.Pupil_plot;
mSM66 = mSM66.Pupil_plot;



%% Before smoothing
HMM_states = [mSM63.HMM_states; mSM64.HMM_states; mSM65.HMM_states; mSM66.HMM_states];
Outcome = [mSM63.Outcome, mSM64.Outcome, mSM65.Outcome, mSM66.Outcome];
Pupil_arousal = [mSM63.Pupil_arousal, mSM64.Pupil_arousal, mSM65.Pupil_arousal, mSM66.Pupil_arousal];

A = Pupil_arousal(Outcome == 0);
B = Pupil_arousal(Outcome == 1);

[~, p_value] = ttest2(A, B);

figure;
subplot(1,2,1);
m = histogram(A, -2:0.05:2, 'DisplayStyle', 'stairs', 'EdgeColor', 'red', 'LineWidth', 1, 'Normalization', 'probability');
hold on
n = histogram(B, -2:0.05:2, 'DisplayStyle', 'stairs', 'EdgeColor', [0 0.4470 0.7410], 'LineWidth', 1, 'Normalization', 'probability');
title(['Dilation(Stimulus - Baseline) vs. Correct/Wrong, p value = ', num2str(p_value)]);
xlim([-2 2]);

y_up = ylim .* 1.1;

line([mean(A) mean(A)], [0 y_up(2)], 'Color', 'red', 'LineWidth', 1, 'LineStyle', '--');
line([mean(B) mean(B)], [0 y_up(2)], 'Color', [0 0.4470 0.7410], 'LineWidth', 1, 'LineStyle', '--');

axis square
xlabel('Pupil Arousal');
ylabel('Probability');

legend([m, n], {'Wrong', 'Correct'});

set(gca,'box','off');
set(gca,'tickdir','out');




[engaged, disengaged] = ArousalState(mSM63.Pupil_arousal, mSM63.HMM_states);

[a, b] = ArousalState(mSM64.Pupil_arousal, mSM64.HMM_states);
engaged = [engaged, a]; disengaged = [disengaged, b];

[a, b] = ArousalState(mSM65.Pupil_arousal, mSM65.HMM_states);
engaged = [engaged, a]; disengaged = [disengaged, b];

[a, b] = ArousalState(mSM66.Pupil_arousal, mSM66.HMM_states);
engaged = [engaged, a]; disengaged = [disengaged, b];



[~, p_value] = ttest2(engaged, disengaged);


subplot(1,2,2);
m = histogram(engaged, -2:0.05:2, 'DisplayStyle', 'stairs', 'EdgeColor', 'red', 'LineWidth', 1, 'Normalization', 'probability');
hold on
n = histogram(disengaged, -2:0.05:2, 'DisplayStyle', 'stairs', 'EdgeColor', [0 0.4470 0.7410], 'LineWidth', 1, 'Normalization', 'probability');
title(['Dilation(Stimulus - Baseline) vs. Engaged/Disengaged, p value = ', num2str(p_value)]);
xlim([-2 2]);

y_up = ylim .* 1.1;

line([mean(engaged) mean(engaged)], [0 y_up(2)], 'Color', 'red', 'LineWidth', 1, 'LineStyle', '--');
line([mean(disengaged) mean(disengaged)], [0 y_up(2)], 'Color', [0 0.4470 0.7410], 'LineWidth', 1, 'LineStyle', '--');

axis square
xlabel('Pupil Arousal');
ylabel('Probability');

legend([m, n], {'Engaged', 'Disengaged'});

set(gca,'box','off');
set(gca,'tickdir','out');




%% After smoothing
mSM63 = SmoothDownsample(mSM63, 50);
mSM64 = SmoothDownsample(mSM64, 50);
mSM65 = SmoothDownsample(mSM65, 50);
mSM66 = SmoothDownsample(mSM66, 50);






HMM_states = [mSM63.HMM_states; mSM64.HMM_states; mSM65.HMM_states; mSM66.HMM_states];
CorrectRate = [mSM63.Outcome, mSM64.Outcome, mSM65.Outcome, mSM66.Outcome];
Pupil_baseline = [mSM63.Pupil_baseline, mSM64.Pupil_baseline, mSM65.Pupil_baseline, mSM66.Pupil_baseline];
Pupil_stim = [mSM63.Pupil_stim, mSM64.Pupil_stim, mSM65.Pupil_stim, mSM66.Pupil_stim];
Pupil_arousal = [mSM63.Pupil_arousal, mSM64.Pupil_arousal, mSM65.Pupil_arousal, mSM66.Pupil_arousal];
% Pupil_arousal = Pupil_stim - Pupil_baseline;




num_sample = length(CorrectRate);
color = [linspace(0,1,num_sample); linspace(0.4470,0,num_sample); linspace(0.7410,0,num_sample)];

[B,I] = sort(HMM_states);





figure;
subplot(2,2,1);
hold on
for i = 1 : num_sample
    
    scatter(CorrectRate(I(i)), Pupil_baseline(I(i)), 30, color(:, i)', 'filled');
    
end

xlim([0.5 1]);
% ylim([-2 2.5]);
% 
% xticks([0.5 1]);
% yticks([-2 0 2.5]);


mdl = fitlm(CorrectRate, Pupil_baseline);
curve1 = plot(xlim, predict(mdl,xlim'));

xlabel('Correct Rate');
ylabel('Pupil Size');
axis square

correlation = corrcoef(CorrectRate,Pupil_baseline);
title(['Baseline epoch, ', num2str(correlation(1,2)), ' / p-value = ', num2str(coefTest(mdl))]);
set(gca,'box','off');
set(gca,'tickdir','out');
hold off





subplot(2,2,2);
hold on
for i = 1 : num_sample
    
    scatter(CorrectRate(I(i)), Pupil_stim(I(i)), 30, color(:, i)', 'filled');
    
end

xlim([0.5 1]);
% ylim([-2 2.5]);`
% 
% xticks([0.5 1]);
% yticks([-2 0 2.5]);


mdl = fitlm(CorrectRate, Pupil_stim);
curve2 = plot(xlim, predict(mdl,xlim'));

xlabel('Correct Rate');
ylabel('Pupil Size');
axis square

correlation = corrcoef(CorrectRate,Pupil_stim);
title(['Stimulus epoch, ', num2str(correlation(1,2)), ' / p-value = ', num2str(coefTest(mdl))]);
set(gca,'box','off');
set(gca,'tickdir','out');
hold off






subplot(2,2,3);
hold on
for i = 1 : num_sample
    
    scatter(CorrectRate(I(i)), Pupil_arousal(I(i)), 30, color(:, i)', 'filled');
    
end

xlim([0.5 1]);


mdl = fitlm(CorrectRate, Pupil_arousal);
curve2 = plot(xlim, predict(mdl,xlim'));

xlabel('Correct Rate');
ylabel('Pupil Arousal');
axis square

correlation = corrcoef(CorrectRate,Pupil_arousal);
title(['r = ', num2str(correlation(1,2)), ' / p-value = ', num2str(coefTest(mdl))]);
set(gca,'box','off');
set(gca,'tickdir','out');
hold off






function data = SmoothDownsample(data, smoothwindow)

data.Pupil_stim = smoothdata(data.Pupil_stim, 'gaussian', smoothwindow);
data.Pupil_baseline = smoothdata(data.Pupil_baseline, 'gaussian', smoothwindow);
data.Pupil_arousal = smoothdata(data.Pupil_arousal, 'gaussian', smoothwindow);
data.Outcome = smoothdata(data.Outcome, 'gaussian', smoothwindow);
data.HMM_states = smoothdata(data.HMM_states, 'gaussian', smoothwindow);


% seed = round(rand * smoothwindow);
seed = 50;


data.Pupil_stim = data.Pupil_stim(seed : smoothwindow : end);
data.Pupil_baseline = data.Pupil_baseline(seed : smoothwindow : end);
data.Pupil_arousal = data.Pupil_arousal(seed : smoothwindow : end);
data.Outcome = data.Outcome(seed : smoothwindow : end);
data.HMM_states = data.HMM_states(seed : smoothwindow : end);

end




function [engaged, disengaged] = ArousalState(arousal, state)

[B,I] = sort(state);

I(isnan(B)) = [];
disengaged_idx = I(1 : floor(length(state) * 0.2));
engaged_idx = I(end - floor(length(state) * 0.2) + 1 : end);


disengaged = arousal(disengaged_idx);
engaged = arousal(engaged_idx);

end
