mSM63 = load('X:\Chaoqun\Papers&Talks\TIV-state Manuscript\Fig5, major_TIM_figure\mSM63_final.mat');
mSM64 = load('X:\Chaoqun\Papers&Talks\TIV-state Manuscript\Fig5, major_TIM_figure\mSM64_final.mat');
mSM65 = load('X:\Chaoqun\Papers&Talks\TIV-state Manuscript\Fig5, major_TIM_figure\mSM65_final.mat');
mSM66 = load('X:\Chaoqun\Papers&Talks\TIV-state Manuscript\Fig5, major_TIM_figure\mSM66_final.mat');
CSP22 = load('X:\Chaoqun\Papers&Talks\TIV-state Manuscript\Fig5, major_TIM_figure\CSP22_final_4state.mat');
CSP23 = load('X:\Chaoqun\Papers&Talks\TIV-state Manuscript\Fig5, major_TIM_figure\CSP23_final_6mice.mat');
CSP32 = load('X:\Chaoqun\Papers&Talks\TIV-state Manuscript\Fig5, major_TIM_figure\CSP32_final_6mice.mat');
CSP38 = load('X:\Chaoqun\Papers&Talks\TIV-state Manuscript\Fig5, major_TIM_figure\CSP38_final.mat');



[all_engaged, all_disengaged] = GroupTrialsEnergy(mSM63);
[e, d] = GroupTrialsEnergy(mSM64);
all_engaged = [all_engaged; e]; all_disengaged = [all_disengaged; d];
[e, d] = GroupTrialsEnergy(mSM65);
all_engaged = [all_engaged; e]; all_disengaged = [all_disengaged; d];
[e, d] = GroupTrialsEnergy(mSM66);
all_engaged = [all_engaged; e]; all_disengaged = [all_disengaged; d];
[e, d] = GroupTrialsEnergy(CSP22);
all_engaged = [all_engaged; e]; all_disengaged = [all_disengaged; d];
[e, d] = GroupTrialsEnergy(CSP23);
all_engaged = [all_engaged; e]; all_disengaged = [all_disengaged; d];
[e, d] = GroupTrialsEnergy(CSP32);
all_engaged = [all_engaged; e]; all_disengaged = [all_disengaged; d];
[e, d] = GroupTrialsEnergy(CSP38);
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


figure('Name', 'Motion energy distribution across states (stimulus+delay), all DLC labels');

% a = histogram(all_engaged, -1:0.02:1.5, 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceAlpha', 0.7);
a = histogram(all_engaged, -1:0.02:1.5, 'DisplayStyle', 'stairs', 'EdgeColor', 'red', 'LineWidth', 1, 'Normalization', 'probability');
hold on
b = histogram(all_disengaged, -1:0.02:1.5, 'DisplayStyle', 'stairs', 'EdgeColor', [0 0.4470 0.7410], 'LineWidth', 1,'Normalization', 'probability');


title({['t test p value = ', num2str(p_value_t)], ['Wilcoxon rank sum test p value = ', num2str(p_value_w)]});

yylim = ylim;
line([mean(all_engaged) mean(all_engaged)], yylim, 'Color', 'red');
line([mean(all_disengaged) mean(all_disengaged)], yylim, 'Color', [0 0.4470 0.7410]);
ylim(yylim);

axis square
xlabel('Normalized Motion Energy');
ylabel('Probability');

legend([a, b], {'Engaged', 'Disengaged'});

set(gca,'box','off');
set(gca,'tickdir','out');

hold off
clear a b









coeff_animal = [mSM63.State_results.StateEnergycorrelation, mSM64.State_results.StateEnergycorrelation,...
    mSM65.State_results.StateEnergycorrelation, mSM66.State_results.StateEnergycorrelation];

coeff_session = [mSM63.State_results.CorreEachSession_motionEnergy, mSM64.State_results.CorreEachSession_motionEnergy,...
    mSM65.State_results.CorreEachSession_motionEnergy, mSM66.State_results.CorreEachSession_motionEnergy];


figure;
hold on

bar([1, 2], [mean(coeff_animal), mean(coeff_session)],'FaceColor',[.8 .8 .8],'EdgeColor',[.3 .3 .3],'LineWidth',1); 
xticks([1, 2]);
set(gca,'xticklabel',{'4 animals', '23 sessions'});
xtickangle(45); 
ylabel('Correlation Coefficient with P(engaged)');

scatter(ones(1,length(coeff_animal)) + (rand(1, length(coeff_animal))-0.5) .* 0.3, coeff_animal, 20, [.5 .5 .5], 'filled');
scatter(ones(1,length(coeff_session)).*2 + (rand(1, length(coeff_session))-0.5) .* 0.3, coeff_session, 20, [.5 .5 .5], 'filled');

er1 = errorbar(1,mean(coeff_animal),std(coeff_animal),std(coeff_animal));    
er1.Color = [0 0 0];                            
er1.LineStyle = 'none';  

er2 = errorbar(2,mean(coeff_session),std(coeff_session),std(coeff_session));    
er2.Color = [0 0 0];                            
er2.LineStyle = 'none';  

set(gca,'box','off');
set(gca,'tickdir','out');
ylim([-1, 1]);
hold off





[mSM63_std_engaged, mSM63_std_disengaged] = MeasureStd(mSM63);
[mSM64_std_engaged, mSM64_std_disengaged] = MeasureStd(mSM64);
[mSM65_std_engaged, mSM65_std_disengaged] = MeasureStd(mSM65);
[mSM66_std_engaged, mSM66_std_disengaged] = MeasureStd(mSM66);



a = [mSM63_std_engaged(1:47); mSM64_std_engaged(1:47); mSM65_std_engaged(1:47); mSM66_std_engaged(1:47)];
b = [mSM63_std_disengaged(1:47); mSM64_std_disengaged(1:47); mSM65_std_disengaged(1:47); mSM66_std_disengaged(1:47)];

figure;
stdshade(a, 0.4, [1 0 0]);
hold on
stdshade(b, 0.4, [0 0.4470 0.7410]);
ylabel('Std of Motion Energy');
ylim([0.2 1]);
set(gca,'box','off');
set(gca,'tickdir','out');



a = [mSM63_std_engaged(end-9 : end); mSM64_std_engaged(end-9 : end); mSM65_std_engaged(end-9 : end); mSM66_std_engaged(end-9 : end)];
b = [mSM63_std_disengaged(end-9 : end); mSM64_std_disengaged(end-9 : end); mSM65_std_disengaged(end-9 : end); mSM66_std_disengaged(end-9 : end)];

figure;
stdshade(a, 0.4, [1 0 0]);
hold on
stdshade(b, 0.4, [0 0.4470 0.7410]);
ylim([0.2 1]);
set(gca,'box','off');
set(gca,'tickdir','out');





function [engaged, disengaged] = GroupTrialsEnergy(data)

% idx_engaged = find(data.State_results.HMM_state > 0.7);
% idx_disengaged = find(data.State_results.HMM_state < 0.3);
% 
% engaged = data.DLCEnergy(idx_engaged);
% disengaged = data.DLCEnergy(idx_disengaged);
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

engaged = data.DLCEnergy(idx_engaged);
disengaged = data.DLCEnergy(idx_disengaged);


end





function [std_engaged, std_disengaged] = MeasureStd(data)

[B,I] = sort(data.State_results.HMM_state);
num_20percent = round(length(I) * 0.2);


idx_engaged = I(end-num_20percent+1 : end);
idx_disengaged = I(1 : num_20percent);




energy_engaged = diff(data.DLC(idx_engaged,:,:), 1, 2);
energy_engaged = sqrt(energy_engaged(:, :, 1:2:end) .^ 2 + energy_engaged(:, :, 2:2:end) .^ 2);

B = permute(energy_engaged,[3 2 1]); 
B_size = size(B);
B = reshape(B, B_size(1), []);

B = normalize(B, 2);
B = reshape(B, B_size);

energy_engaged = permute(B,[3 2 1]);
energy_engaged = nanmean(energy_engaged, 3);
std_engaged = nanstd(energy_engaged, 1);
clear B B_size


energy_disengaged = diff(data.DLC(idx_disengaged,:,:), 1, 2);
energy_disengaged = sqrt(energy_disengaged(:, :, 1:2:end) .^ 2 + energy_disengaged(:, :, 2:2:end) .^ 2);

B = permute(energy_disengaged,[3 2 1]); 
B_size = size(B);
B = reshape(B, B_size(1), []);

B = normalize(B, 2);
B = reshape(B, B_size);

energy_disengaged = permute(B,[3 2 1]);
energy_disengaged = nanmean(energy_disengaged, 3);
std_disengaged = nanstd(energy_disengaged, 1);
clear B B_size


end
