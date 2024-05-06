mSM63 = load('X:\Chaoqun\Papers&Talks\TIV-state Manuscript\TIM_without_outcomeRegressor_Fig\mSM63_NoOutcomeRegressor.mat');
mSM64 = load('X:\Chaoqun\Papers&Talks\TIV-state Manuscript\TIM_without_outcomeRegressor_Fig\mSM64_NoOutcomeRegressor.mat');
mSM65 = load('X:\Chaoqun\Papers&Talks\TIV-state Manuscript\TIM_without_outcomeRegressor_Fig\mSM65_NoOutcomeRegressor.mat');
mSM66 = load('X:\Chaoqun\Papers&Talks\TIV-state Manuscript\TIM_without_outcomeRegressor_Fig\mSM66_NoOutcomeRegressor.mat');



coeff_TIM = [mSM63.State_results.StateTIMcorrelation, mSM64.State_results.StateTIMcorrelation, ...
    mSM65.State_results.StateTIMcorrelation, mSM66.State_results.StateTIMcorrelation];

coeff_TIM_sessions = [mSM63.State_results.CorreEachSession_TIM, mSM64.State_results.CorreEachSession_TIM, ...
    mSM65.State_results.CorreEachSession_TIM, mSM66.State_results.CorreEachSession_TIM];


c = corrcoef(mSM63.DLCEnergy, mSM63.CorrectRate');
d = corrcoef(mSM64.DLCEnergy, mSM64.CorrectRate');
e = corrcoef(mSM65.DLCEnergy, mSM65.CorrectRate');
f = corrcoef(mSM66.DLCEnergy, mSM66.CorrectRate');


coeff_motionEnergy = [c(1,2),d(1,2),e(1,2),f(1,2)];



%%
figure('name', 'Bar graph');
hold on

bar([1, 2], [mean(coeff_TIM), nanmean(coeff_TIM_sessions)],'FaceColor',[.8 .8 .8],'EdgeColor',[.3 .3 .3],'LineWidth',1); 
xticks([1, 2]);
set(gca,'xticklabel',{'4 animals', '23 sessions'});
xtickangle(45); 
ylabel('Correlation Coefficient');

scatter(ones(1,length(coeff_TIM)) + (rand(1, length(coeff_TIM))-0.5) .* 0.3, coeff_TIM, 20, [0 0.4470 0.7410], 'filled');
scatter(ones(1,length(coeff_TIM_sessions)).*2 + (rand(1, length(coeff_TIM_sessions))-0.5) .* 0.3, coeff_TIM_sessions, 20, [0 0.4470 0.7410], 'filled');

er1 = errorbar(1,mean(coeff_TIM),std(coeff_TIM),std(coeff_TIM));    
er1.Color = [0 0 0];                            
er1.LineStyle = 'none';  

er2 = errorbar(2,nanmean(coeff_TIM_sessions),nanstd(coeff_TIM_sessions),nanstd(coeff_TIM_sessions));    
er2.Color = [0 0 0];                            
er2.LineStyle = 'none';  

set(gca,'box','off');
set(gca,'tickdir','out');
ylim([-0.8, 0.4]);
hold off



%%
a = mSM66.Distance_TIV(mSM66.borders(4) : mSM66.borders(5));
b = mSM66.State_results.HMM_state(mSM66.borders(4) : mSM66.borders(5));


figure('Name', 'Example sessions, mSM66');

j = plot(a);
hold on
ylabel('Standardized TIM');

yyaxis right
k = plot(b);
ylabel('P(engaged)');
ylim([0 1]);
xlabel('Trial Number');
xlim([0 length(a)]);



[R, P] = corrcoef(a, b);
title(['Correcoef = ', num2str(R(1,2)), ', p = ', num2str(P(1,2))]);
legend([j, k], {'TIM', 'P(engaged)'});
set(gca,'box','off');
set(gca,'tickdir','out');

xticks([0 length(a)]);

hold off




%%
TIM = [mSM63.Distance_TIV; mSM64.Distance_TIV; mSM65.Distance_TIV; mSM66.Distance_TIV];
correctRate = [mSM63.CorrectRate, mSM64.CorrectRate, mSM65.CorrectRate, mSM66.CorrectRate];
states = [mSM63.State_results.HMM_state; mSM64.State_results.HMM_state; mSM65.State_results.HMM_state; mSM66.State_results.HMM_state];



a = TIM(34:50:end);
b = correctRate(34:50:end);
c = states(34:50:end);

[B,I] = sort(c);

num_sample = length(a);
color = [linspace(0,1,num_sample); linspace(0.4470,0,num_sample); linspace(0.7410,0,num_sample)];



figure;
hold on
for i = 1 : num_sample
    
    scatter(b(I(i)), a(I(i)), 20, color(:, i)', 'filled');
    
end

xlim([0.5 1]);
ylim([-0.3 0.5]);

Fit_linear = polyfit(b, a', 1);
curve1 = plot(xlim, polyval(Fit_linear,xlim));

[R,P] = corrcoef(b, a');


xlabel('Correct Rate');
ylabel('TIM');
title(['CorrCoef = ', num2str(R(1,2)), ', p = ', num2str(P(1,2))]);
axis square

set(gca,'box','off');
set(gca,'tickdir','out');
hold off
