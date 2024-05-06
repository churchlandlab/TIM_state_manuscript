CSP38 = load('X:\Chaoqun\Papers&Talks\TIV-state Manuscript\Fig5, major_TIM_figure\CSP38_final.mat');
CSP32 = load('X:\Chaoqun\Papers&Talks\TIV-state Manuscript\Fig5, major_TIM_figure\CSP32_final.mat');
CSP23 = load('X:\Chaoqun\Papers&Talks\TIV-state Manuscript\Fig5, major_TIM_figure\CSP23_final.mat');
CSP22 = load('X:\Chaoqun\Papers&Talks\TIV-state Manuscript\Fig5, major_TIM_figure\CSP22_final.mat');

Fez13 = load('X:\Chaoqun\Papers&Talks\TIV-state Manuscript\Fig5, major_TIM_figure\Fez13_final.mat');
Fez7 = load('X:\Chaoqun\Papers&Talks\TIV-state Manuscript\Fig5, major_TIM_figure\Fez7_final.mat');
CTP7 = load('X:\Chaoqun\Papers&Talks\TIV-state Manuscript\Fig5, major_TIM_figure\CTP7_final.mat');
Plex05 = load('X:\Chaoqun\Papers&Talks\TIV-state Manuscript\Fig5, major_TIM_figure\Plex05_final.mat');

mSM85 = load('X:\Chaoqun\Papers&Talks\TIV-state Manuscript\Fig5, major_TIM_figure\mSM85_final.mat');




coeff_animal = [CSP22.Results.TIVPerformanceCoef, CSP23.Results.TIVPerformanceCoef, CSP32.Results.TIVPerformanceCoef, CSP38.Results.TIVPerformanceCoef,...
    Fez13.Results.TIVPerformanceCoef, Fez7.Results.TIVPerformanceCoef, CTP7.Results.TIVPerformanceCoef, Plex05.Results.TIVPerformanceCoef, ...
    mSM85.Results.TIVPerformanceCoef];

coeff_session = [CSP22.Results.CorreEachSession_TIV, CSP23.Results.CorreEachSession_TIV, CSP32.Results.CorreEachSession_TIV, CSP38.Results.CorreEachSession_TIV,...
    Fez13.Results.CorreEachSession_TIV, Fez7.Results.CorreEachSession_TIV, CTP7.Results.CorreEachSession_TIV, Plex05.Results.CorreEachSession_TIV, ...
    mSM85.Results.CorreEachSession_TIV];





subplot(1,2,1);
hold on

bar([1, 2], [mean(coeff_animal), mean(coeff_session)],'FaceColor',[.8 .8 .8],'EdgeColor',[.3 .3 .3],'LineWidth',1); 
xticks([1, 2]);
set(gca,'xticklabel',{'9 animals', '73 sessions'});
xtickangle(45); 
ylabel('Correlation Coefficient with Correct Rate');

scatter(ones(1,length(coeff_animal)) + (rand(1, length(coeff_animal))-0.5) .* 0.4, coeff_animal, 20, [0 0.4470 0.7410], 'filled');
scatter(ones(1,length(coeff_session)).*2 + (rand(1, length(coeff_session))-0.5) .* 0.4, coeff_session, 20, [0 0.4470 0.7410], 'filled');

er1 = errorbar(1,mean(coeff_animal),std(coeff_animal),std(coeff_animal));    
er1.Color = [0 0 0];                            
er1.LineStyle = 'none';  

er2 = errorbar(2,mean(coeff_session),std(coeff_session),std(coeff_session));    
er2.Color = [0 0 0];                            
er2.LineStyle = 'none';  

set(gca,'box','off');
set(gca,'tickdir','out');
ylim([-1, 0.4]);
hold off






coeff_animal = [CSP22.Results.EnergyPerformanceCoef(401), CSP23.Results.EnergyPerformanceCoef(401), CSP32.Results.EnergyPerformanceCoef(401), CSP38.Results.EnergyPerformanceCoef(401),...
    Fez13.Results.EnergyPerformanceCoef(201), Fez7.Results.EnergyPerformanceCoef(401), CTP7.Results.EnergyPerformanceCoef(201), Plex05.Results.EnergyPerformanceCoef(201), ...
    mSM85.Results.EnergyPerformanceCoef(201)];

coeff_session = [CSP22.Results.CorreEachSession_motionEnergy, CSP23.Results.CorreEachSession_motionEnergy, CSP32.Results.CorreEachSession_motionEnergy, CSP38.Results.CorreEachSession_motionEnergy,...
    Fez13.Results.CorreEachSession_motionEnergy, Fez7.Results.CorreEachSession_motionEnergy, CTP7.Results.CorreEachSession_motionEnergy, Plex05.Results.CorreEachSession_motionEnergy, ...
    mSM85.Results.CorreEachSession_motionEnergy];


subplot(1,2,2);
hold on

bar([1, 2], [mean(coeff_animal), mean(coeff_session)],'FaceColor',[.8 .8 .8],'EdgeColor',[.3 .3 .3],'LineWidth',1); 
xticks([1, 2]);
set(gca,'xticklabel',{'9 animals', '73 sessions'});
xtickangle(45); 
ylabel('Correlation Coefficient with Correct Rate');

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
