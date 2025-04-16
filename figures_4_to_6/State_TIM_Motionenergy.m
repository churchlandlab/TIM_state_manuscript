function [Results] = State_TIM_Motionenergy(HMM_state, TIM, motionEnergy, session_borders, smoothwindow, raw_TIM)

global mouse_name
mousename = mouse_name;

%% 1st, Plot HMM state (smoothed), TIM, and motion energy
figure('Name', ['P(engaged), TIM, motion energy, ', mousename]);
hold on
a = plot(TIM);
b = plot(motionEnergy);

ylimit = [max([TIM; motionEnergy; HMM_state]), min([TIM; motionEnergy; HMM_state])] .* 1.2;
for i = 1 : (length(session_borders) - 1)
    line([session_borders(i) session_borders(i)], ylimit, 'LineStyle','--', 'Color', [.5 .5 .5]);
end
ylabel('TIM/Motion Energy');

yyaxis right
c = plot(HMM_state);
ylabel('Probability of Engaged State');
xlabel('Trial Number');
legend([a, b, c], {'TIM', 'Motion Energy', 'P(engaged)'});
set(gca,'box','off');
set(gca,'tickdir','out');
hold off
clear a b c




%% 2nd, HMM state-TIM correlation per session
coeff_session = nan(1, length(session_borders));
coeff_motionEnergy = nan(1, length(session_borders));

ttt = [1, session_borders];

for i = 1 : length(session_borders)
    a = HMM_state(ttt(i): ttt(i+1)-1);
    b = TIM(ttt(i): ttt(i+1)-1);
    c = motionEnergy(ttt(i): ttt(i+1)-1);
    
    x = corrcoef(a, b);
    coeff_session(i) = x(1,2);
    
    xx = corrcoef(a, c);
    coeff_motionEnergy(i) = xx(1,2);
end


Results.CorreEachSession_TIM = coeff_session;
Results.CorreEachSession_motionEnergy = coeff_motionEnergy;



figure('Name', ['Single-session correlation coefficient, ', mousename]);
hold on

bar([1, 2], [mean(coeff_session), mean(coeff_motionEnergy)],'FaceColor',[.7 .7 .7],'EdgeColor',[.3 .3 .3],'LineWidth',1); 
xticks([1, 2]);
set(gca,'xticklabel',{'TIM', 'Motion Energy'});
xtickangle(45); 
ylabel('Correlation Coefficient');

scatter(ones(1,length(coeff_session)) + (rand(1, length(coeff_session))-0.5) .* 0.3, coeff_session, 20, [0 0.4470 0.7410], 'filled');
scatter(ones(1,length(coeff_session)).*2 + (rand(1, length(coeff_session))-0.5) .* 0.3, coeff_motionEnergy, 20, [.3 .3 .3], 'filled');

er1 = errorbar(1,mean(coeff_session),std(coeff_session),std(coeff_session));    
er1.Color = [0 0 0];                            
er1.LineStyle = 'none';  

er2 = errorbar(2,mean(coeff_motionEnergy),std(coeff_motionEnergy),std(coeff_motionEnergy));    
er2.Color = [0 0 0];                            
er2.LineStyle = 'none';  

set(gca,'box','off');
set(gca,'tickdir','out');
hold off

clear i a b c x xx er1 er2 ttt





%% 3rd, Linear shift test (+-200 trials)
corre_shift = nan(1, 401);
corre_shift_2 = nan(1, 401);
for i = -200 : 200
    TIM_shifted = circshift(TIM, i);
    motionEnergy_shifted = circshift(motionEnergy, i);
    a = corrcoef(TIM_shifted, HMM_state);
    b = corrcoef(motionEnergy_shifted, HMM_state);
    corre_shift(i+201) = a(1,2);  
    corre_shift_2(i+201) = b(1,2); 
end


Results.CorreShift = corre_shift;
Results.EnergyPerformanceCoef = corre_shift_2;

figure('Name', ['P(engaged)-TIM/motion Energy correlation temporal shift, ', mousename]);
hold on
plot(corre_shift);
plot(corre_shift_2);
xlim([1 401]);
xticks([1 101 201 301 401]);
xticklabels({'-200','-100','0','100','200'});
line([201 201], ylim, 'Color','black','LineStyle','--');
xlabel('Trial Shift');
ylabel('Correlation Coefficient');
legend({'TIM', 'Motion Energy'});
hold off
clear i a b






%% 4th, Trend of performance and TIM changes within each session
% We normalize all sessions with more than 300 trials to 500-trial long.
ttt = [1, session_borders];
HMM_state_sessions = [];
TIM_sessions = [];

for i = 1 : length(session_borders)
    a = HMM_state(ttt(i): ttt(i+1)-1);
    b = TIM(ttt(i): ttt(i+1)-1);
    
    if length(a) > 300 && length(a) == length(b)
        zoom = 1:500/length(a):500;
        if length(zoom) < length(a)   % sometimes you get an array 1 element shorter than array a.
            zoom = [zoom, 500];
        end
        
        a_adjusted = interp1(zoom, a, 1:500, 'linear');    
        b_adjusted = interp1(zoom, b, 1:500, 'linear');
        
        HMM_state_sessions = [HMM_state_sessions; a_adjusted];
        TIM_sessions = [TIM_sessions; b_adjusted];
        
    end
    
end


HMM_state_sessions(:, end) = [];
figure('Name', ['TIM/P(state) fluctuations within each session, ', mousename]);
subplot(2, 1, 1);
hold on
for i = 1 : size(HMM_state_sessions, 1)
    plot(HMM_state_sessions(i, :), 'Color', [0.4 0.4 0.4]);
end
xticks([1 250 500]);
xticklabels({'Session Beginning','Session Middle','Session End'});
ylabel('P(engaged)');

autocorrelation = (sum(sum(corrcoef(HMM_state_sessions'))) - size(HMM_state_sessions,1)) / ...
    (size(HMM_state_sessions,1)*(size(HMM_state_sessions,1)-1));
title(['Cross-session correlation coeff: ', num2str(autocorrelation)]);





TIM_sessions(:, end) = [];
subplot(2, 1, 2);
hold on
for i = 1 : size(TIM_sessions, 1)
    plot(TIM_sessions(i, :), 'Color', [0.4 0.4 0.4]);
end
xticks([1 250 500]);
xticklabels({'Session Beginning','Session Middle','Session End'});
ylabel('TIM');

autocorrelation = (sum(sum(corrcoef(TIM_sessions'))) - size(TIM_sessions,1)) / ...
    (size(TIM_sessions,1)*(size(TIM_sessions,1)-1));
title(['Cross-session correlation coeff: ', num2str(autocorrelation)]);

clear a b zoom a_adjusted b_adjusted ttt i



x = corrcoef(HMM_state_sessions');
y = corrcoef(TIM_sessions');

for i = 1 : size(x, 1)
    x(i,i) = NaN;
    y(i,i) = NaN;
end
x = reshape(x, 1, []);
y = reshape(y, 1, []);

Results.CrossSessionCorrelation_state = x;
Results.CrossSessionCorrelation_TIM = y;


end