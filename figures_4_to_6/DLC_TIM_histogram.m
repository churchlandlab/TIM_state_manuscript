mSM63 = load('C:\Users\churchland\Desktop\State-TIM paper double check\mSM63_shared.mat');
mSM64 = load('C:\Users\churchland\Desktop\State-TIM paper double check\mSM64_shared.mat');
mSM65 = load('C:\Users\churchland\Desktop\State-TIM paper double check\mSM65_shared.mat');
mSM66 = load('C:\Users\churchland\Desktop\State-TIM paper double check\mSM66_shared.mat');




[all_engaged, all_disengaged] = GroupTrialsTIM(mSM63);
animal_label = ones(length(all_engaged), 1);

[e, d] = GroupTrialsTIM(mSM64);
all_engaged = [all_engaged; e]; all_disengaged = [all_disengaged; d];
animal_label = [animal_label; ones(length(e), 1).*2];

[e, d] = GroupTrialsTIM(mSM65);
all_engaged = [all_engaged; e]; all_disengaged = [all_disengaged; d];
animal_label = [animal_label; ones(length(e), 1).*3];

[e, d] = GroupTrialsTIM(mSM66);
all_engaged = [all_engaged; e]; all_disengaged = [all_disengaged; d];
animal_label = [animal_label; ones(length(e), 1).*4];





%% Mixed-effects model test
group1 = all_engaged; 
group2 = all_disengaged; 

% Combine data into a single table
data = [group1; group2]; % Combined data
group = [repmat({'Group1'}, size(group1), 1); ...
         repmat({'Group2'}, size(group2), 1)]; % Group labels
subjectID = [animal_label; animal_label]; % Subject IDs

subjectID_1 = subjectID == 1; 
subjectID_2 = subjectID == 2;
subjectID_3 = subjectID == 3;
subjectID_4 = subjectID == 4;


% Create a table for the mixed-effects model
tbl = table(data, group, subjectID_1, subjectID_2, subjectID_3, subjectID_4,...
    'VariableNames', {'Data', 'Group', 'SubjectID_1', 'SubjectID_2', 'SubjectID_3', 'SubjectID_4'});

% Define the mixed-effects model
% Random intercept for each subject (SubjectID)
% Fixed effect for Group
% lme = fitlme(tbl, 'Data ~ Group + (1|SubjectID)');
lme = fitlme(tbl, 'Data ~ Group + (1|SubjectID_1) + (1|SubjectID_2) + (1|SubjectID_3) + (1|SubjectID_4)');

% Display results
disp(lme);






%% Histogram
figure('Name', 'TIM distribution across states (stimulus+delay), all DLC labels');

a = histogram(all_engaged, -1:0.02:1.5, 'DisplayStyle', 'stairs', 'EdgeColor', 'red', 'LineWidth', 1, 'Normalization', 'probability');
hold on
b = histogram(all_disengaged, -1:0.02:1.5, 'DisplayStyle', 'stairs', 'EdgeColor', [0 0.4470 0.7410], 'LineWidth', 1,'Normalization', 'probability');

yylim = ylim;
line([mean(all_engaged) mean(all_engaged)], yylim, 'Color', 'red');
line([mean(all_disengaged) mean(all_disengaged)], yylim, 'Color', [0 0.4470 0.7410]);
ylim(yylim);

axis square
xlabel('Normalized TIM');
ylabel('Probability');

legend([a, b], {'Engaged', 'Disengaged'});

set(gca,'box','off');
set(gca,'tickdir','out');

hold off
clear a b




function [engaged, disengaged] = GroupTrialsTIM(data)

[B,I] = sort(data.State_results.HMM_state);
num_20percent = round(length(I) * 0.2);


idx_engaged = I(end-num_20percent+1 : end);
idx_disengaged = I(1 : num_20percent);

engaged = data.TIM(idx_engaged);
disengaged = data.TIM(idx_disengaged);


end
