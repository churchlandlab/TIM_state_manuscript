function Pupil_plot = Pupil_performance_state(fillPupil, raw_data, HMM_states, smoothwindow)

CorrectRate = raw_data.Rewarded;
CorrectRate = smoothdata(CorrectRate, 'gaussian', smoothwindow);



Pupil_baseline = cell2mat(fillPupil.Baseline);
Pupil_stim = cell2mat(fillPupil.StimuWindow);

Pupil_stim = mean(Pupil_stim, 1, 'omitnan'); 
Pupil_baseline = mean(Pupil_baseline, 1, 'omitnan'); 
Pupil_arousal = Pupil_stim - Pupil_baseline;



Pupil_plot.Pupil_stim = Pupil_stim;
Pupil_plot.Pupil_baseline = Pupil_baseline;
Pupil_plot.Pupil_arousal = Pupil_arousal;
Pupil_plot.Outcome = raw_data.Rewarded;
Pupil_plot.HMM_states = HMM_states;




Pupil_stim = smoothdata(Pupil_stim, 'gaussian', smoothwindow);
Pupil_baseline = smoothdata(Pupil_baseline, 'gaussian', smoothwindow);

seed = round(rand * smoothwindow);
CorrectRate = CorrectRate(seed : smoothwindow : end);
Pupil_stim = Pupil_stim(seed : smoothwindow : end);
Pupil_baseline = Pupil_baseline(seed : smoothwindow : end);
HMM_states = HMM_states(seed : smoothwindow : end);



%% Plot correct rate vs. pupil size, with color code denoting P(engaged)
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
% ylim([-2 2.5]);
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

A = Pupil_arousal(raw_data.Rewarded == 0);
B = Pupil_arousal(raw_data.Rewarded == 1);

[~, p_value] = ttest2(A, B);


m = histogram(A, -2:0.05:3, 'DisplayStyle', 'stairs', 'EdgeColor', 'red', 'LineWidth', 1, 'Normalization', 'probability');
hold on
n = histogram(B, -2:0.05:3, 'DisplayStyle', 'stairs', 'EdgeColor', [0 0.4470 0.7410], 'LineWidth', 1, 'Normalization', 'probability');
title(['Dilation(Stimulus - Baseline) vs. Correct/Wrong, p value = ', num2str(p_value)]);
xlim([-2 3]);

y_up = ylim .* 1.1;

line([mean(A) mean(A)], [0 y_up(2)], 'Color', 'red', 'LineWidth', 1, 'LineStyle', '--');
line([mean(B) mean(B)], [0 y_up(2)], 'Color', [0 0.4470 0.7410], 'LineWidth', 1, 'LineStyle', '--');

axis square
xlabel('Pupil Arousal');
ylabel('Probability');

legend([m, n], {'Wrong', 'Correct'});

set(gca,'box','off');
set(gca,'tickdir','out');





end