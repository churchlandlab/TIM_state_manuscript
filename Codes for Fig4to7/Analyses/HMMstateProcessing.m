%% This function is designed to match the trial numbers of GLM-HMM state file and processed DLC/behaivoral files
%% This issue is generated because some trials are removed when calculating the state and preparing the DLC/behavioral files
%% These steps are done separtely. We use this function to put them on the same page.

function engaged_state = HMMstateProcessing(statefile, DLC_deleteIdx, raw_data)

HMM_deleteIdx = cell2mat(statefile.masks(1:length(statefile.masks)));
HMM_deleteIdx = find(HMM_deleteIdx~=1);

HMM_prop = cell2mat(statefile.posterior_probs(1:length(statefile.masks))');


% There are 3 columns in the probs varaible. We need to know which one denotes the engaged state
% There is a "state_label_indices" variable in the HMM file, it contains 3 numbers, 
% they denote the column number of engaged/left bias/right bias states in the the probs varaible
idx_engaged = statefile.state_label_indices(1);
idx_leftbias = statefile.state_label_indices(2);
idx_rightbias = statefile.state_label_indices(3);


for i = 1 : length(HMM_deleteIdx)
    a = HMM_prop(1 : HMM_deleteIdx(i)-1, :);
    b = HMM_prop(HMM_deleteIdx(i):end, :);
    
    HMM_prop = cat(1, a, nan(1,3), b);
end


HMM_prop(DLC_deleteIdx, :) = [];

engaged_state = HMM_prop(:, idx_engaged);
leftbias_state = HMM_prop(:, idx_leftbias);
rightbias_state = HMM_prop(:, idx_rightbias);


trials_engaged = find(engaged_state >= 0.8);
trials_leftbias = find(leftbias_state >= 0.8);
trials_rightbias = find(rightbias_state >= 0.8);


figure('Name', 'The choice pattern in different HMM states');
subplot(1,3,1);
hold on
a = raw_data.CorrectSide(trials_engaged);
b = raw_data.ResponseSide(trials_engaged);

plot([1,2], [sum(a == 1 & a == b)/sum(a == 1), sum(a == 2 & a == b)/sum(a == 2)], 'Color', [0 0.4470 0.7410]);
xlim([0.7 2.3]);
xticks([1, 2]);
xticklabels({'Left','Right'});
ylim([0, 1]);

xlabel('Correct Side');
ylabel('Correct Rate');
title('Engaged');
set(gca,'box','off');

set(gca,'tickdir','out');
hold off



subplot(1,3,2);
hold on
a = raw_data.CorrectSide(trials_leftbias);
b = raw_data.ResponseSide(trials_leftbias);

plot([1,2], [sum(a == 1 & a == b)/sum(a == 1), sum(a == 2 & a == b)/sum(a == 2)], 'Color', [0.8500 0.3250 0.0980]);
xlim([0.7 2.3]);
xticks([1, 2]);
xticklabels({'Left','Right'});
ylim([0, 1]);

title('Left-biased');
set(gca,'box','off');

set(gca,'tickdir','out');
hold off


subplot(1,3,3);
hold on
a = raw_data.CorrectSide(trials_rightbias);
b = raw_data.ResponseSide(trials_rightbias);

plot([1,2], [sum(a == 1 & a == b)/sum(a == 1), sum(a == 2 & a == b)/sum(a == 2)], 'Color', [0.9290 0.6940 0.1250]);
xlim([0.7 2.3]);
xticks([1, 2]);
xticklabels({'Left','Right'});
ylim([0, 1]);

title('Right-biased');
set(gca,'box','off');

set(gca,'tickdir','out');
hold off


end