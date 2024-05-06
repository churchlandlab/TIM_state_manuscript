SessionList(9,:)
SessionList(29,:)
SessionList(42,:)
SessionList(58,:)
SessionList(71,:)
SessionList(92,:)

mSM64_state = cell2mat(posterior_probs');

plot(mSM64_state(:,1));
xlabel('Trial Num');
ylabel('P(engaged)');



SessionList([50,51,53,54,56],:)




masks = masks([50,51,53,54,56]);
model_training_sessions = model_training_sessions([50,51,53,54,56],:);
mouse = mouse([50,51,53,54,56],:);
posterior_probs = posterior_probs([50,51,53,54,56]);


subplot(2,3,6);
hold on
scatter(mSM63.State_results.CrossSessionCorrelation_state, mSM63.State_results.CrossSessionCorrelation_TIM, 40, 'filled');
scatter(mSM64.State_results.CrossSessionCorrelation_state, mSM64.State_results.CrossSessionCorrelation_TIM, 40, 'filled');
scatter(mSM65.State_results.CrossSessionCorrelation_state, mSM65.State_results.CrossSessionCorrelation_TIM, 40, 'filled');
scatter(mSM66.State_results.CrossSessionCorrelation_state, mSM66.State_results.CrossSessionCorrelation_TIM, 40, 'filled');

xlim([-1 1]);
ylim([-0.4 1]);
line([0 0], ylim, 'Color', 'k', 'LineStyle', '--');
line(xlim, [0 0], 'Color', 'k', 'LineStyle', '--');

xlabel('Cross-session Corr, State');
ylabel('Cross-session Corr, TIM');
set(gca,'box','off');
set(gca,'tickdir','out');
axis square
hold off









coeff_animal = [mSM63.State_results.StateTIMcorrelation, mSM64.State_results.StateTIMcorrelation,...
    mSM65.State_results.StateTIMcorrelation, mSM66.State_results.StateTIMcorrelation];

coeff_session = [mSM63.State_results.CorreEachSession_TIM, mSM64.State_results.CorreEachSession_TIM,...
    mSM65.State_results.CorreEachSession_TIM, mSM66.State_results.CorreEachSession_TIM];


subplot(2,3,2);
hold on

bar([1, 2], [mean(coeff_animal), mean(coeff_session)],'FaceColor',[.8 .8 .8],'EdgeColor',[.3 .3 .3],'LineWidth',1); 
xticks([1, 2]);
set(gca,'xticklabel',{'4 animals', '23 sessions'});
xtickangle(45); 
ylabel('Correlation Coefficient with P(engaged)');

scatter(ones(1,length(coeff_animal)) + (rand(1, length(coeff_animal))-0.5) .* 0.3, coeff_animal, 20, [0 0.4470 0.7410], 'filled');
scatter(ones(1,length(coeff_session)).*2 + (rand(1, length(coeff_session))-0.5) .* 0.3, coeff_session, 20, [0 0.4470 0.7410], 'filled');

er1 = errorbar(1,mean(coeff_animal),std(coeff_animal),std(coeff_animal));    
er1.Color = [0 0 0];                            
er1.LineStyle = 'none';  

er2 = errorbar(2,mean(coeff_session),std(coeff_session),std(coeff_session));    
er2.Color = [0 0 0];                            
er2.LineStyle = 'none';  

set(gca,'box','off');
set(gca,'tickdir','out');
ylim([-0.8, 0.4]);
hold off







subplot(2,3,5);
hold on
plot(mSM63.State_results.CorreShift);
plot(mSM64.State_results.CorreShift);
plot(mSM65.State_results.CorreShift);
plot(mSM66.State_results.CorreShift);

ylim([-0.5 0.3]);

xlim([0 400]);
line([200 200], ylim, 'Color', 'k', 'LineStyle', '--');
xticks([0, 200, 400]);
xticklabels({'-200','0','200'});
xlabel('Trial Shift');
ylabel('Correlation Coefficient ');
set(gca,'box','off');
set(gca,'tickdir','out');

axis square
hold off






mSM63_sigma = (mSM63.State_results.StateTIMcorrelation - mean(mSM63.State_results.CoefSessionShuffle)) / std(mSM63.State_results.CoefSessionShuffle);
mSM64_sigma = (mSM64.State_results.StateTIMcorrelation - mean(mSM64.State_results.CoefSessionShuffle)) / std(mSM64.State_results.CoefSessionShuffle);
mSM65_sigma = (mSM65.State_results.StateTIMcorrelation - mean(mSM65.State_results.CoefSessionShuffle)) / std(mSM65.State_results.CoefSessionShuffle);
mSM66_sigma = (mSM66.State_results.StateTIMcorrelation - mean(mSM66.State_results.CoefSessionShuffle)) / std(mSM66.State_results.CoefSessionShuffle);

control_mean = mean(mSM63.State_results.CoefSessionShuffle);
control_std = std(mSM63.State_results.CoefSessionShuffle);
control_Distribution = mSM63.State_results.CoefSessionShuffle ./ control_std;


subplot(2,3,4);
hold on
h = histogram(control_Distribution, 'FaceColor', [0 0.4470 0.7410], 'EdgeColor', 'none', 'FaceAlpha', 0.8);
line([mSM63_sigma mSM63_sigma], ylim, 'Color','red');
line([mSM64_sigma mSM64_sigma], ylim, 'Color','red');
line([mSM65_sigma mSM65_sigma], ylim, 'Color','red');
line([mSM66_sigma mSM66_sigma], ylim, 'Color','red');


x = min(control_Distribution) : 0.001 : max(control_Distribution);
y = normpdf(x,mean(control_Distribution),std(control_Distribution));

ttt = max(h.Values) / max(y);
plot(x,y.*ttt,'Color', [0.6 0.6 0.6], 'LineWidth', 2);


xlabel('Correlation Coefficient (Standard Deviation)');
ylabel('Number of Shuffles');

set(gca,'box','off');
set(gca,'tickdir','out');

axis square
hold off





CorrectRate = [smoothdata(mSM63.Outcome, 'gaussian', 50), smoothdata(mSM64.Outcome, 'gaussian', 50), ...
    smoothdata(mSM65.Outcome, 'gaussian', 50), smoothdata(mSM66.Outcome, 'gaussian', 50)];
TIM = [smoothdata(mSM63.TIM, 'gaussian', 50); smoothdata(mSM64.TIM, 'gaussian', 50); ...
    smoothdata(mSM65.TIM, 'gaussian', 50); smoothdata(mSM66.TIM, 'gaussian', 50)];
HMM = [mSM63.State_results.HMM_state; mSM64.State_results.HMM_state; mSM65.State_results.HMM_state; mSM66.State_results.HMM_state];


CorrectRate = CorrectRate(42:50:end);
TIM = TIM(42:50:end);
HMM = HMM(42:50:end);



num_sample = length(CorrectRate);
color = [linspace(0,1,num_sample); linspace(0.4470,0,num_sample); linspace(0.7410,0,num_sample)];

[B,I] = sort(HMM);


subplot(2,3,3);
hold on
for i = 1 : num_sample
    
    scatter(CorrectRate(I(i)), TIM(I(i)), 30, color(:, i)', 'filled');
    
end

xlim([0.5 1]);
ylim([-0.4 0.8]);
% 
% xticks([0.5 1]);
% yticks([-2 0 2.5]);


mdl = fitlm(CorrectRate, TIM);
curve2 = plot(xlim, predict(mdl,xlim'));

xlabel('Correct Rate');
ylabel('Standardized TIM');
axis square

correlation = corrcoef(CorrectRate,TIM);
title(['Corr = ', num2str(correlation(1,2)), ' / p-value = ', num2str(coefTest(mdl))]);
set(gca,'box','off');
set(gca,'tickdir','out');
axis square
hold off






a = [mSM63.State_results.CrossSessionCorrelation_state, mSM64.State_results.CrossSessionCorrelation_state,...
    mSM65.State_results.CrossSessionCorrelation_state, mSM66.State_results.CrossSessionCorrelation_state];

b = [mSM63.State_results.CrossSessionCorrelation_TIM, mSM64.State_results.CrossSessionCorrelation_TIM,...
    mSM65.State_results.CrossSessionCorrelation_TIM, mSM66.State_results.CrossSessionCorrelation_TIM];


figure;
histogram(a, -1:0.1:1, 'Orientation', 'vertical', 'FaceColor', [.5 .5 .5]);
line([nanmean(a),nanmean(a)], ylim, 'LineStyle','--', 'Color', [1 0 0]);
xlim([-1 1])
set(gca,'box','off');
set(gca,'tickdir','out');


figure;
histogram(b, -0.4:0.1:1, 'Orientation', 'horizontal', 'FaceColor', [.5 .5 .5]);
line(xlim, [nanmean(b),nanmean(b)], 'LineStyle','--', 'Color', [1 0 0]);
ylim([-0.4 1])
set(gca,'box','off');
set(gca,'tickdir','out');



