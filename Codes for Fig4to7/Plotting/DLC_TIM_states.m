mSM63 = load('X:\Chaoqun\Papers&Talks\TIV-state Manuscript\Fig5, major_TIM_figure\mSM63_final.mat');
mSM64 = load('X:\Chaoqun\Papers&Talks\TIV-state Manuscript\Fig5, major_TIM_figure\mSM64_final.mat');
mSM65 = load('X:\Chaoqun\Papers&Talks\TIV-state Manuscript\Fig5, major_TIM_figure\mSM65_final.mat');
mSM66 = load('X:\Chaoqun\Papers&Talks\TIV-state Manuscript\Fig5, major_TIM_figure\mSM66_final.mat');
CSP22 = load('X:\Chaoqun\Papers&Talks\TIV-state Manuscript\Fig5, major_TIM_figure\CSP22_final.mat');
CSP38 = load('X:\Chaoqun\Papers&Talks\TIV-state Manuscript\Fig5, major_TIM_figure\CSP38_final.mat');



[all_engaged, all_disengaged] = GroupTrialsTIM(mSM63);
[e, d] = GroupTrialsTIM(mSM64);
all_engaged = [all_engaged; e]; all_disengaged = [all_disengaged; d];
[e, d] = GroupTrialsTIM(mSM65);
all_engaged = [all_engaged; e]; all_disengaged = [all_disengaged; d];
[e, d] = GroupTrialsTIM(mSM66);
all_engaged = [all_engaged; e]; all_disengaged = [all_disengaged; d];
[e, d] = GroupTrialsTIM(CSP22);
all_engaged = [all_engaged; e]; all_disengaged = [all_disengaged; d];
[e, d] = GroupTrialsTIM(CSP23);
all_engaged = [all_engaged; e]; all_disengaged = [all_disengaged; d];
[e, d] = GroupTrialsTIM(CSP32);
all_engaged = [all_engaged; e]; all_disengaged = [all_disengaged; d];
[e, d] = GroupTrialsTIM(CSP38);
all_engaged = [all_engaged; e]; all_disengaged = [all_disengaged; d];


% if length(all_disengaged) > length(all_engaged)
%     
%     idx = randperm(length(all_disengaged), length(all_engaged)); 
%     all_disengaged = all_disengaged(idx);
%     
% else
%     idx = randperm(length(all_engaged), length(all_disengaged));
%     all_engaged = all_engaged(idx);
%     
% end



[~, p_value_t] = ttest2(all_engaged, all_disengaged);
[p_value_w, ~] = ranksum(all_engaged, all_disengaged);


figure('Name', 'TIM distribution across states (stimulus+delay), all DLC labels');

a = histogram(all_engaged, -1:0.02:1.5, 'DisplayStyle', 'stairs', 'EdgeColor', 'red', 'LineWidth', 1, 'Normalization', 'probability');
hold on
b = histogram(all_disengaged, -1:0.02:1.5, 'DisplayStyle', 'stairs', 'EdgeColor', [0 0.4470 0.7410], 'LineWidth', 1,'Normalization', 'probability');

title({['t test p value = ', num2str(p_value_t)], ['Wilcoxon rank sum test p value = ', num2str(p_value_w)]});

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

% idx_engaged = find(data.State_results.HMM_state > 0.7);
% idx_disengaged = find(data.State_results.HMM_state < 0.3);
% 
% engaged = data.TIM(idx_engaged);
% disengaged = data.TIM(idx_disengaged);
% 
% if length(disengaged) > length(engaged)
%     
%     idx = randperm(length(disengaged), length(engaged)); 
%     disengaged = disengaged(idx);
%     
% else
%     idx = randperm(length(engaged), length(disengaged));
%     engaged = engaged(idx);
%     
% end

[B,I] = sort(data.State_results.HMM_state);
num_20percent = round(length(I) * 0.2);


idx_engaged = I(end-num_20percent+1 : end);
idx_disengaged = I(1 : num_20percent);

engaged = data.TIM(idx_engaged);
disengaged = data.TIM(idx_disengaged);


end
