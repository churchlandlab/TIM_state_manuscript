load('X:\Chaoqun\Papers&Talks\TIV-state Manuscript\Pupil_fig\JC047.mat');
load('X:\Chaoqun\Papers&Talks\TIV-state Manuscript\Pupil_fig\JC066.mat');
load('X:\Chaoqun\Papers&Talks\TIV-state Manuscript\Pupil_fig\JC072.mat');
load('X:\Chaoqun\Papers&Talks\TIV-state Manuscript\Pupil_fig\JC099.mat');
load('X:\Chaoqun\Papers&Talks\TIV-state Manuscript\Pupil_fig\JC103.mat');
load('X:\Chaoqun\Papers&Talks\TIV-state Manuscript\Pupil_fig\JC111.mat');



%%
a = corrcoef(JC047.Pupil, JC047.Performance, 'Rows','complete');
b = corrcoef(JC066.Pupil, JC066.Performance, 'Rows','complete');
c = corrcoef(JC072.Pupil, JC072.Performance, 'Rows','complete');
d = corrcoef(JC099.Pupil, JC099.Performance, 'Rows','complete');
e = corrcoef(JC103.Pupil, JC103.Performance, 'Rows','complete');
f = corrcoef(JC111.Pupil, JC111.Performance, 'Rows','complete');

corr_pupil = [a(1,2), b(1,2), c(1,2), d(1,2), e(1,2), f(1,2)];



a = corrcoef(JC047.TIM, JC047.Performance, 'Rows','complete');
b = corrcoef(JC066.TIM, JC066.Performance, 'Rows','complete');
c = corrcoef(JC072.TIM, JC072.Performance, 'Rows','complete');
d = corrcoef(JC099.TIM, JC099.Performance, 'Rows','complete');
e = corrcoef(JC103.TIM, JC103.Performance, 'Rows','complete');
f = corrcoef(JC111.TIM, JC111.Performance, 'Rows','complete');

corr_TIM = [a(1,2), b(1,2), c(1,2), d(1,2), e(1,2), f(1,2)];


a = corrcoef(JC047.DLCEnergy, JC047.Performance, 'Rows','complete');
b = corrcoef(JC066.DLCEnergy, JC066.Performance, 'Rows','complete');
c = corrcoef(JC072.DLCEnergy, JC072.Performance, 'Rows','complete');
d = corrcoef(JC099.DLCEnergy, JC099.Performance, 'Rows','complete');
e = corrcoef(JC103.DLCEnergy, JC103.Performance, 'Rows','complete');
f = corrcoef(JC111.DLCEnergy, JC111.Performance, 'Rows','complete');

corr_Energy = [a(1,2), b(1,2), c(1,2), d(1,2), e(1,2), f(1,2)];


figure('name', '6 animals');
hold on

bar([1, 2, 3], [mean(corr_pupil), mean(corr_Energy), mean(corr_TIM)],'FaceColor',[.8 .8 .8],'EdgeColor',[.3 .3 .3],'LineWidth',1); 
xticks([1, 2, 3]);
set(gca,'xticklabel',{'Pupil', 'Motion Energy', 'TIM'});
xtickangle(45); 
ylabel('Correlation Coefficient');

scatter(ones(1,length(corr_pupil)) + (rand(1, length(corr_pupil))-0.5) .* 0.3, corr_pupil, 20, [.5 .5 .5], 'filled');
scatter(ones(1,length(corr_Energy)).*2 + (rand(1, length(corr_Energy))-0.5) .* 0.3, corr_Energy, 20, [.5 .5 .5], 'filled');
scatter(ones(1,length(corr_TIM)).*3 + (rand(1, length(corr_TIM))-0.5) .* 0.3, corr_TIM, 20, [.5 .5 .5], 'filled');

er1 = errorbar(1,mean(corr_pupil),std(corr_pupil),std(corr_pupil));    
er1.Color = [0 0 0];                            
er1.LineStyle = 'none';  


er2 = errorbar(2,mean(corr_Energy),std(corr_Energy),std(corr_Energy));    
er2.Color = [0 0 0];                            
er2.LineStyle = 'none';  


er3 = errorbar(3,mean(corr_TIM),std(corr_TIM),std(corr_TIM));    
er3.Color = [0 0 0];                            
er3.LineStyle = 'none';  

set(gca,'box','off');
set(gca,'tickdir','out');
hold off





%%
a = normalize(JC047.Pupil);
a = a(25:50:end);
b = normalize(JC066.Pupil);
b = b(25:50:end);
c = normalize(JC072.Pupil);
c = c(25:50:end);
d = normalize(JC099.Pupil);
d = d(25:50:end);
e = normalize(JC103.Pupil);
e = e(25:50:end);
f = normalize(JC111.Pupil);
f = f(25:50:end);

figure('Name', 'Pupil vs. Performance, Joao Dataset');
hold on
scatter(JC047.Performance(25:50:end), a, 20, 'filled');
scatter(JC066.Performance(25:50:end), b, 20, 'filled');
scatter(JC072.Performance(25:50:end), c, 20, 'filled');
scatter(JC099.Performance(25:50:end), d, 20, 'filled');
scatter(JC103.Performance(25:50:end), e, 20, 'filled');
scatter(JC111.Performance(25:50:end), f, 20, 'filled');

axis square
xlim([0.5 1]);


m = [JC047.Performance(25:50:end); JC066.Performance(25:50:end); JC072.Performance(25:50:end); ...
    JC099.Performance(25:50:end); JC103.Performance(25:50:end); JC111.Performance(25:50:end)];
n = [a; b; c; d; e; f];

Fit_linear = polyfit(m, n, 1);
curve1 = plot(xlim, polyval(Fit_linear,xlim));

[coef, pvalue] = corrcoef(m, n);

xlabel('Correct Rate');
ylabel('Standardized Pupil Size');

title(['CorrCoef = ', num2str(coef(1,2)), '; p-value = ', num2str(pvalue(1, 2))]);
set(gca,'box','off');
set(gca,'tickdir','out');
hold off







