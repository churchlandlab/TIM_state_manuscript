mSM63 = load('C:\Users\churchland\Desktop\State-TIM paper double check\mSM63_results.mat');
mSM64 = load('C:\Users\churchland\Desktop\State-TIM paper double check\mSM64_results.mat');
mSM65 = load('C:\Users\churchland\Desktop\State-TIM paper double check\mSM65_results.mat');
mSM66 = load('C:\Users\churchland\Desktop\State-TIM paper double check\mSM66_results.mat');


%% Mixed-effects Model
%% TIM
a = mSM63.State_results.CorreEachSession_TIM;
b = mSM64.State_results.CorreEachSession_TIM;
c = mSM65.State_results.CorreEachSession_TIM;
d = mSM66.State_results.CorreEachSession_TIM;
% Example data
animal = [ones(1, length(a)), 2 * ones(1, length(b)), 3 * ones(1, length(c)), 4 * ones(1, length(d))]'; % Animal IDs
y = [a, b, c, d]';


% Combine into a table
data = table(animal, y, 'VariableNames', {'Animals', 'TIM'});
disp(data);

% Define the formula for the mixed-effects model
formula = 'TIM ~ 1 + (1|Animals)'; % Test if the overall mean is different from 0, accounting for group variability
% Fit the model
lme = fitlme(data, formula);
% Display the results
disp(lme);



%% Motion Energy
a = mSM63.State_results.CorreEachSession_motionEnergy;
b = mSM64.State_results.CorreEachSession_motionEnergy;
c = mSM65.State_results.CorreEachSession_motionEnergy;
d = mSM66.State_results.CorreEachSession_motionEnergy;
% Example data
y = [a, b, c, d]';


% Combine into a table
data = table(animal, y, 'VariableNames', {'Animals', 'Motion'});
disp(data);

% Define the formula for the mixed-effects model
formula = 'Motion ~ 1 + (1|Animals)'; % Test if the overall mean is different from 0, accounting for group variability
% Fit the model
lme = fitlme(data, formula);
% Display the results
disp(lme);



%% Cross-session correlation, P(engaged)
a = mSM63.State_results.CrossSessionCorrelation_state;
b = mSM64.State_results.CrossSessionCorrelation_state;
c = mSM65.State_results.CrossSessionCorrelation_state;
d = mSM66.State_results.CrossSessionCorrelation_state;

animal = [ones(1, sum(~isnan(a))), 2 * ones(1, sum(~isnan(b))), 3 * ones(1, sum(~isnan(c))), 4 * ones(1, sum(~isnan(d)))]'; % Animal IDs
y = [a, b, c, d]';
y(isnan(y)) = [];

% Combine into a table
data = table(animal, y, 'VariableNames', {'Animals', 'Cross_Session_State'});
disp(data);

% Define the formula for the mixed-effects model
formula = 'Cross_Session_State ~ 1 + (1|Animals)'; % Test if the overall mean is different from 0, accounting for group variability
% Fit the model
lme = fitlme(data, formula);
% Display the results
disp(lme);



%% Cross-session correlation TIM
a = mSM63.State_results.CrossSessionCorrelation_TIM;
b = mSM64.State_results.CrossSessionCorrelation_TIM;
c = mSM65.State_results.CrossSessionCorrelation_TIM;
d = mSM66.State_results.CrossSessionCorrelation_TIM;

animal = [ones(1, sum(~isnan(a))), 2 * ones(1, sum(~isnan(b))), 3 * ones(1, sum(~isnan(c))), 4 * ones(1, sum(~isnan(d)))]'; % Animal IDs
y = [a, b, c, d]';
y(isnan(y)) = [];

% Combine into a table
data = table(animal, y, 'VariableNames', {'Animals', 'Cross_Session_State'});
disp(data);

% Define the formula for the mixed-effects model
formula = 'Cross_Session_State ~ 1 + (1|Animals)'; % Test if the overall mean is different from 0, accounting for group variability
% Fit the model
lme = fitlme(data, formula);
% Display the results
disp(lme);












%% Bar Plot
%% TIM
a = mSM63.State_results.CorreEachSession_TIM;
b = mSM64.State_results.CorreEachSession_TIM;
c = mSM65.State_results.CorreEachSession_TIM;
d = mSM66.State_results.CorreEachSession_TIM;


figure('Name', 'Correlation coeff, all sessions');
subplot(1, 2, 1);
hold on

scatter(ones(1,length(a)) + (rand(1, length(a))-0.5) .* 0.2, a, 40, [0.4660 0.6740 0.1880], 'filled');
scatter(ones(1,length(b)) + (rand(1, length(b))-0.5) .* 0.2, b, 40, [0.8500 0.3250 0.0980], 'filled');
scatter(ones(1,length(c)) + (rand(1, length(c))-0.5) .* 0.2, c, 40, [0.9290 0.6940 0.1250], 'filled');
scatter(ones(1,length(d)) + (rand(1, length(d))-0.5) .* 0.2, d, 40, [0.4940 0.1840 0.5560], 'filled');
xticks([1]);
xticklabels({'TIM'});

er1 = errorbar(1,mean([a,b,c,d]),std([a,b,c,d]),std([a,b,c,d]));    
er1.Color = [0 0 0];                            
er1.LineStyle = 'none';  

line([0.8 1.2], [mean([a,b,c,d]) mean([a,b,c,d])], 'Color', [0 0 0], 'linewidth',2);
line([0.5 1.5], [0 0], 'Color', [0.5 0.5 0.5], 'linewidth', 0.5, 'LineStyle', '--');

ylabel('Correlation Coefficient with P(engaged)');

set(gca,'box','off'); set(gca,'tickdir','out');
xlim([0.5 1.5]);
ylim([-1, 1]);
hold off






%% Motion Energy
a = mSM63.State_results.CorreEachSession_motionEnergy;
b = mSM64.State_results.CorreEachSession_motionEnergy;
c = mSM65.State_results.CorreEachSession_motionEnergy;
d = mSM66.State_results.CorreEachSession_motionEnergy;


subplot(1, 2, 2);
hold on

scatter(ones(1,length(a)) + (rand(1, length(a))-0.5) .* 0.2, a, 40, [0.4660 0.6740 0.1880], 'filled');
scatter(ones(1,length(b)) + (rand(1, length(b))-0.5) .* 0.2, b, 40, [0.8500 0.3250 0.0980], 'filled');
scatter(ones(1,length(c)) + (rand(1, length(c))-0.5) .* 0.2, c, 40, [0.9290 0.6940 0.1250], 'filled');
scatter(ones(1,length(d)) + (rand(1, length(d))-0.5) .* 0.2, d, 40, [0.4940 0.1840 0.5560], 'filled');
xticks([1]);
xticklabels({'Motion Energy'});

er1 = errorbar(1,mean([a,b,c,d]),std([a,b,c,d]),std([a,b,c,d]));    
er1.Color = [0 0 0];                            
er1.LineStyle = 'none';  

line([0.8 1.2], [mean([a,b,c,d]) mean([a,b,c,d])], 'Color', [0 0 0], 'linewidth',2);
line([0.5 1.5], [0 0], 'Color', [0.5 0.5 0.5], 'linewidth', 0.5, 'LineStyle', '--');

ylabel('Correlation Coefficient with P(engaged)');

set(gca,'box','off'); set(gca,'tickdir','out');
xlim([0.5 1.5]);
ylim([-1, 1]);
hold off




%% Cross-session correlations
figure;
hold on
scatter(mSM63.State_results.CrossSessionCorrelation_state, mSM63.State_results.CrossSessionCorrelation_TIM, 40, [0.4660 0.6740 0.1880], 'filled');
scatter(mSM64.State_results.CrossSessionCorrelation_state, mSM64.State_results.CrossSessionCorrelation_TIM, 40, [0.8500 0.3250 0.0980], 'filled');
scatter(mSM65.State_results.CrossSessionCorrelation_state, mSM65.State_results.CrossSessionCorrelation_TIM, 40, [0.9290 0.6940 0.1250], 'filled');
scatter(mSM66.State_results.CrossSessionCorrelation_state, mSM66.State_results.CrossSessionCorrelation_TIM, 40, [0.4940 0.1840 0.5560], 'filled');

xlim([-1 1]);
ylim([0 1]);
line([0 0], ylim, 'Color', 'k', 'LineStyle', '--');

xlabel('Cross-session Corr, State');
ylabel('Cross-session Corr, TIM');
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
histogram(b, 0:0.05:1, 'Orientation', 'horizontal', 'FaceColor', [.5 .5 .5]);
line(xlim, [nanmean(b),nanmean(b)], 'LineStyle','--', 'Color', [1 0 0]);
ylim([0 1])
set(gca,'box','off');
set(gca,'tickdir','out');





