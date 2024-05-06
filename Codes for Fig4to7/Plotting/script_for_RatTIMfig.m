load('X:\Chaoqun\Papers&Talks\TIV-state Manuscript\Rat TIM fig\cy11.mat');
load('X:\Chaoqun\Papers&Talks\TIV-state Manuscript\Rat TIM fig\cy39.mat');
load('X:\Chaoqun\Papers&Talks\TIV-state Manuscript\Rat TIM fig\cy42.mat');
load('X:\Chaoqun\Papers&Talks\TIV-state Manuscript\Rat TIM fig\cy44.mat');
load('X:\Chaoqun\Papers&Talks\TIV-state Manuscript\Rat TIM fig\cy47.mat');


%%
a = cy44.TIM;
b = cy44.Performance;

a = smoothdata(a, 'gaussian', 50);
b = smoothdata(b, 'gaussian', 50);

figure('Name', 'Example sessions, cy44');

j = plot(a);
hold on
ylabel('Standardized TIM');

yyaxis right
k = plot(b);
ylabel('Correct Rate');
ylim([0.3 1]);
xlabel('Trial Number');

for i = 2 : (length(cy44.SessionBorders) - 1)
    line([cy44.SessionBorders(i) cy44.SessionBorders(i)], [0.3 1], 'Color', 'black', 'LineStyle', ':');
end

correlation = corrcoef(a, b);
title(['Correcoef = ', num2str(correlation(1,2))]);
legend([j, k], {'TIV', 'Performance'});
set(gca,'box','off');
set(gca,'tickdir','out');
hold off






%%
coeff_TIM = [cy11.CorrCoeff_TIM, cy39.CorrCoeff_TIM, cy42.CorrCoeff_TIM, cy44.CorrCoeff_TIM, cy47.CorrCoeff_TIM];
coeff_motionEnergy = [cy11.CorrCoeff_motion, cy39.CorrCoeff_motion, cy42.CorrCoeff_motion, cy44.CorrCoeff_motion, cy47.CorrCoeff_motion];



coeff_TIM_sessions = [];
coeff_motionEnergy_sessions = [];

for i = 1 : 2
    a = cy11.TIM(cy11.SessionBorders(i) + 1 : cy11.SessionBorders(i+1));
    b = cy11.Performance(cy11.SessionBorders(i) + 1 : cy11.SessionBorders(i+1));
    c = cy11.motionEnergy(cy11.SessionBorders(i) + 1 : cy11.SessionBorders(i+1));
    
    a = smoothdata(a, 'gaussian', 50);
    b = smoothdata(b, 'gaussian', 50); 
    c = smoothdata(c, 'gaussian', 50);
    
    u = corrcoef(a, b);
    v = corrcoef(c, b);
    
    coeff_TIM_sessions = [coeff_TIM_sessions, u(1,2)];
    coeff_motionEnergy_sessions = [coeff_motionEnergy_sessions, v(1,2)];
    
end



for i = 1 : 2
    a = cy39.TIM(cy39.SessionBorders(i) + 1 : cy39.SessionBorders(i+1));
    b = cy39.Performance(cy39.SessionBorders(i) + 1 : cy39.SessionBorders(i+1));
    c = cy39.motionEnergy(cy39.SessionBorders(i) + 1 : cy39.SessionBorders(i+1));
    
    a = smoothdata(a, 'gaussian', 50);
    b = smoothdata(b, 'gaussian', 50); 
    c = smoothdata(c, 'gaussian', 50);
    
    u = corrcoef(a, b);
    v = corrcoef(c, b);
    
    coeff_TIM_sessions = [coeff_TIM_sessions, u(1,2)];
    coeff_motionEnergy_sessions = [coeff_motionEnergy_sessions, v(1,2)];
    
end


for i = 1 : 3
    a = cy42.TIM(cy42.SessionBorders(i) + 1 : cy42.SessionBorders(i+1));
    b = cy42.Performance(cy42.SessionBorders(i) + 1 : cy42.SessionBorders(i+1));
    c = cy42.motionEnergy(cy42.SessionBorders(i) + 1 : cy42.SessionBorders(i+1));
    
    a = smoothdata(a, 'gaussian', 50);
    b = smoothdata(b, 'gaussian', 50); 
    c = smoothdata(c, 'gaussian', 50);
    
    u = corrcoef(a, b);
    v = corrcoef(c, b);
    
    coeff_TIM_sessions = [coeff_TIM_sessions, u(1,2)];
    coeff_motionEnergy_sessions = [coeff_motionEnergy_sessions, v(1,2)];
    
end


for i = 1 : 3
    a = cy44.TIM(cy44.SessionBorders(i) + 1 : cy44.SessionBorders(i+1));
    b = cy44.Performance(cy44.SessionBorders(i) + 1 : cy44.SessionBorders(i+1));
    c = cy44.motionEnergy(cy44.SessionBorders(i) + 1 : cy44.SessionBorders(i+1));
    
    a = smoothdata(a, 'gaussian', 50);
    b = smoothdata(b, 'gaussian', 50); 
    c = smoothdata(c, 'gaussian', 50);
    
    u = corrcoef(a, b);
    v = corrcoef(c, b);
    
    coeff_TIM_sessions = [coeff_TIM_sessions, u(1,2)];
    coeff_motionEnergy_sessions = [coeff_motionEnergy_sessions, v(1,2)];
    
end


for i = 1 : 4
    a = cy47.TIM(cy47.SessionBorders(i) + 1 : cy47.SessionBorders(i+1));
    b = cy47.Performance(cy47.SessionBorders(i) + 1 : cy47.SessionBorders(i+1));
    c = cy47.motionEnergy(cy47.SessionBorders(i) + 1 : cy47.SessionBorders(i+1));
    
    a = smoothdata(a, 'gaussian', 50);
    b = smoothdata(b, 'gaussian', 50); 
    c = smoothdata(c, 'gaussian', 50);
    
    u = corrcoef(a, b);
    v = corrcoef(c, b);
    
    coeff_TIM_sessions = [coeff_TIM_sessions, u(1,2)];
    coeff_motionEnergy_sessions = [coeff_motionEnergy_sessions, v(1,2)];
    
end





figure('name', 'All sessions');
subplot(1,2,1);
hold on

bar([1, 2], [mean(coeff_TIM), mean(coeff_motionEnergy)],'FaceColor',[.8 .8 .8],'EdgeColor',[.3 .3 .3],'LineWidth',1); 
xticks([1, 2]);
set(gca,'xticklabel',{'TIM', 'Motion Energy'});
xtickangle(45); 
ylabel('Correlation Coefficient');

scatter(ones(1,length(coeff_TIM)) + (rand(1, length(coeff_TIM))-0.5) .* 0.3, coeff_TIM, 20, [0 0.4470 0.7410], 'filled');
scatter(ones(1,length(coeff_TIM)).*2 + (rand(1, length(coeff_TIM))-0.5) .* 0.3, coeff_motionEnergy, 20, [.5 .5 .5], 'filled');

er1 = errorbar(1,mean(coeff_TIM),std(coeff_TIM),std(coeff_TIM));    
er1.Color = [0 0 0];                            
er1.LineStyle = 'none';  

er2 = errorbar(2,mean(coeff_motionEnergy),std(coeff_motionEnergy),std(coeff_motionEnergy));    
er2.Color = [0 0 0];                            
er2.LineStyle = 'none';  

set(gca,'box','off');
set(gca,'tickdir','out');
ylim([-1, 0.6]);
hold off



subplot(1,2,2);
hold on

bar([1, 2], [mean(coeff_TIM_sessions), mean(coeff_motionEnergy_sessions)],'FaceColor',[.8 .8 .8],'EdgeColor',[.3 .3 .3],'LineWidth',1); 
xticks([1, 2]);
set(gca,'xticklabel',{'TIM', 'Motion Energy'});
xtickangle(45); 
ylabel('Correlation Coefficient');

scatter(ones(1,length(coeff_TIM_sessions)) + (rand(1, length(coeff_TIM_sessions))-0.5) .* 0.3, coeff_TIM_sessions, 20, [0 0.4470 0.7410], 'filled');
scatter(ones(1,length(coeff_TIM_sessions)).*2 + (rand(1, length(coeff_TIM_sessions))-0.5) .* 0.3, coeff_motionEnergy_sessions, 20, [.5 .5 .5], 'filled');

er1 = errorbar(1,mean(coeff_TIM_sessions),std(coeff_TIM_sessions),std(coeff_TIM_sessions));    
er1.Color = [0 0 0];                            
er1.LineStyle = 'none';  

er2 = errorbar(2,mean(coeff_motionEnergy_sessions),std(coeff_motionEnergy_sessions),std(coeff_motionEnergy_sessions));    
er2.Color = [0 0 0];                            
er2.LineStyle = 'none';  

set(gca,'box','off');
set(gca,'tickdir','out');
ylim([-1, 0.6]);
hold off





%%
figure('Name', 'Shift test, rats');
hold on

plot(cy11.TIV_shift);
plot(cy39.TIV_shift);
plot(cy42.TIV_shift);
plot(cy44.TIV_shift);
plot(cy47.TIV_shift);

xlim([1 401]);
ylim([-0.8 0.4]);

line([200 200], [-0.8 0.4], 'Color', 'black', 'LineStyle', ':');

xlabel('Trial Shift');
ylabel('Correlation Coefficient');

set(gca,'box','off');
set(gca,'tickdir','out');
hold off







