mSM63 = load('X:\Chaoqun\Papers&Talks\TIV-state Manuscript\Fig5, major_TIM_figure\mSM63_final.mat');
mSM64 = load('X:\Chaoqun\Papers&Talks\TIV-state Manuscript\Fig5, major_TIM_figure\mSM64_final.mat');
mSM65 = load('X:\Chaoqun\Papers&Talks\TIV-state Manuscript\Fig5, major_TIM_figure\mSM65_final.mat');
mSM66 = load('X:\Chaoqun\Papers&Talks\TIV-state Manuscript\Fig5, major_TIM_figure\mSM66_final.mat');
% mSM85 = load('X:\Chaoqun\Papers&Talks\TIV-state Manuscript\Fig5, major_TIM_figure\mSM85_final.mat');
% Fez7 = load('X:\Chaoqun\Papers&Talks\TIV-state Manuscript\Fig5, major_TIM_figure\Fez7_final.mat');
% Fez13 = load('X:\Chaoqun\Papers&Talks\TIV-state Manuscript\Fig5, major_TIM_figure\Fez13_final.mat');
% CTP7 = load('X:\Chaoqun\Papers&Talks\TIV-state Manuscript\Fig5, major_TIM_figure\CTP7_final.mat');
% Plex05 = load('X:\Chaoqun\Papers&Talks\TIV-state Manuscript\Fig5, major_TIM_figure\Plex05_final.mat');
% CSP22 = load('X:\Chaoqun\Papers&Talks\TIV-state Manuscript\Fig5, major_TIM_figure\CSP22_final.mat');
% CSP38 = load('X:\Chaoqun\Papers&Talks\TIV-state Manuscript\Fig5, major_TIM_figure\CSP38_final.mat');


A = []; B = []; C = []; D = [];
[a,b,c,d] = EachAnimal(mSM63.Outcome, mSM63.TIM, mSM63.State_results.HMM_state);
A = [A; a]; B = [B; b]; C = [C; c]; D = [D; d];
[a,b,c,d] = EachAnimal(mSM64.Outcome, mSM64.TIM, mSM64.State_results.HMM_state);
A = [A; a]; B = [B; b]; C = [C; c]; D = [D; d];
[a,b,c,d] = EachAnimal(mSM65.Outcome, mSM65.TIM, mSM65.State_results.HMM_state);
A = [A; a]; B = [B; b]; C = [C; c]; D = [D; d];
[a,b,c,d] = EachAnimal(mSM66.Outcome, mSM66.TIM, mSM66.State_results.HMM_state);
A = [A; a]; B = [B; b]; C = [C; c]; D = [D; d];

% [a,b,c,d] = EachAnimal(mSM85.Outcome, mSM85.TIM);
% A = [A; a]; B = [B; b]; C = [C; c]; D = [D; d];
% 
% [a,b,c,d] = EachAnimal(Fez7.Outcome, Fez7.TIM);
% A = [A; a]; B = [B; b]; C = [C; c]; D = [D; d];
% 
% [a,b,c,d] = EachAnimal(Fez13.Outcome, Fez13.TIM);
% A = [A; a]; B = [B; b]; C = [C; c]; D = [D; d];
% 
% [a,b,c,d] = EachAnimal(CTP7.Outcome, CTP7.TIM);
% A = [A; a]; B = [B; b]; C = [C; c]; D = [D; d];
% 
% [a,b,c,d] = EachAnimal(Plex05.Outcome, Plex05.TIM);
% A = [A; a]; B = [B; b]; C = [C; c]; D = [D; d];
% 
% [a,b,c,d] = EachAnimal(CSP22.Outcome, CSP22.TIM);
% A = [A; a]; B = [B; b]; C = [C; c]; D = [D; d];
% 
% [a,b,c,d] = EachAnimal(CSP38.Outcome, CSP38.TIM);
% A = [A; a]; B = [B; b]; C = [C; c]; D = [D; d];



[~, p_value] = ttest2(A,D);
figure;
subplot(1,2,1);


m = histogram(A, -1:0.05:1.3, 'DisplayStyle', 'stairs', 'EdgeColor', 'red', 'LineWidth', 1, 'Normalization', 'probability');
hold on
n = histogram(D, -1:0.05:1.3, 'DisplayStyle', 'stairs', 'EdgeColor', [0 0.4470 0.7410], 'LineWidth', 1,'Normalization', 'probability');
title(['Outcome(n-1) = 0, p value = ', num2str(p_value)]);
xlim([-1 1.3]);

y_up = ylim .* 1.1;

line([mean(A) mean(A)], [0 y_up(2)], 'Color', 'red', 'LineWidth', 1, 'LineStyle', '--');
line([mean(D) mean(D)], [0 y_up(2)], 'Color', [0 0.4470 0.7410], 'LineWidth', 1, 'LineStyle', '--');

axis square
xlabel('TIM');
ylabel('Number of Trials');

legend([m, n], {'Disengaged', 'Engaged'});

set(gca,'box','off');
set(gca,'tickdir','out');

hold off
clear m n





[~, p_value] = ttest2(C, B);
subplot(1,2,2);
m = histogram(C, -1:0.05:1.3, 'DisplayStyle', 'stairs', 'EdgeColor', 'red', 'LineWidth', 1, 'Normalization', 'probability');
hold on
n = histogram(B, -1:0.05:1.3, 'DisplayStyle', 'stairs', 'EdgeColor', [0 0.4470 0.7410], 'LineWidth', 1,'Normalization', 'probability');
title(['Outcome(n-1) = 1, p value = ', num2str(p_value)]);
xlim([-1 1.3]);

y_up = ylim .* 1.1;

line([mean(c) mean(c)], [0 y_up(2)], 'Color', 'red', 'LineWidth', 1, 'LineStyle', '--');
line([mean(b) mean(b)], [0 y_up(2)], 'Color', [0 0.4470 0.7410], 'LineWidth', 1, 'LineStyle', '--');

axis square
xlabel('TIM');
ylabel('Number of Trials');

legend([m, n], {'Disengaged', 'Engaged'});

set(gca,'box','off');
set(gca,'tickdir','out');

hold off
clear m n





function [a,b,c,d] = EachAnimal(Outcome, TIM_raw, HMM)

Current = Outcome;
Previous = [0, Outcome];
Previous(end) = [];

bar_engaged = sort(HMM);
bar_disengaged = bar_engaged(round(length(bar_engaged) * 0.2));
bar_engaged = bar_engaged(round(length(bar_engaged) * 0.8));



% a = TIM_raw(Current == 0 & Previous == 0);
% b = TIM_raw(Current == 1 & Previous == 1);
% c = TIM_raw(Current == 0 & Previous == 1);
% d = TIM_raw(Current == 1 & Previous == 0);

a = TIM_raw(HMM <= bar_disengaged & Previous' == 0);
b = TIM_raw(HMM >= bar_engaged & Previous' == 1);
c = TIM_raw(HMM <= bar_disengaged & Previous' == 1);
d = TIM_raw(HMM >= bar_engaged & Previous' == 0);



if length(a) > length(d)
    p = randperm(length(a), length(d));
    a = a(p);
    
elseif length(a) < length(d)
    p = randperm(length(d), length(a));
    d = d(p);
end


if length(b) > length(c)
    p = randperm(length(b), length(c));
    b = b(p);
    
elseif length(b) < length(c)
    p = randperm(length(c), length(b));
    c = c(p);
end



end
