CTP7_x = CTP7.Results.TIVPerformanceCoef / std(CTP7.Results.CoefSessionShuffle);
Plex05_x = CTP7.Results.TIVPerformanceCoef / std(Plex05.Results.CoefSessionShuffle);
Fez7_x = Fez7.Results.TIVPerformanceCoef / std(Fez7.Results.CoefSessionShuffle);
Fez13_x = Fez13.Results.TIVPerformanceCoef / std(Fez13.Results.CoefSessionShuffle);
mSM63_x = mSM63.Results.TIVPerformanceCoef / std(mSM63.Results.CoefSessionShuffle);
mSM64_x = mSM64.Results.TIVPerformanceCoef / std(mSM64.Results.CoefSessionShuffle);
mSM65_x = mSM65.Results.TIVPerformanceCoef / std(mSM65.Results.CoefSessionShuffle);
mSM66_x = mSM66.Results.TIVPerformanceCoef / std(mSM66.Results.CoefSessionShuffle);
mSM85_x = mSM85.Results.TIVPerformanceCoef / std(mSM85.Results.CoefSessionShuffle);


plot_std = std(Fez7.Results.CoefSessionShuffle);
plot_mean = mean(Fez7.Results.CoefSessionShuffle);
coef_shuffle = Fez7.Results.CoefSessionShuffle;

figure;
hold on
m = histogram(Fez7.Results.CoefSessionShuffle - plot_mean, 'FaceColor', [0 0.4470 0.7410], 'EdgeColor', 'none', 'FaceAlpha', 0.8);
a = line([CTP7_x CTP7_x] .* plot_std, ylim, 'Color','red');
b = line([Plex05_x Plex05_x] .* plot_std, ylim, 'Color','red');
c = line([Fez7_x Fez7_x] .* plot_std, ylim, 'Color','red');
d = line([Fez13_x Fez13_x] .* plot_std, ylim, 'Color','red');
e = line([mSM63_x mSM63_x] .* plot_std, ylim, 'Color','red');
f = line([mSM64_x mSM64_x] .* plot_std, ylim, 'Color','red');
g = line([mSM65_x mSM65_x] .* plot_std, ylim, 'Color','red');
h = line([mSM66_x mSM66_x] .* plot_std, ylim, 'Color','red');
i = line([mSM85_x mSM85_x] .* plot_std, ylim, 'Color','red');



x = min(coef_shuffle) : 0.001 : max(coef_shuffle);
y = normpdf(x,mean(coef_shuffle),std(coef_shuffle));

ttt = max(m.Values) / max(y);
plot(x,y.*ttt,'Color', [0.6 0.6 0.6], 'LineWidth', 2);

legend(a, {'True Correlation Coeff'});
xlabel('Distance to Shuffle Control (sigma)');
ylabel('Number of Shuffles');


xticks([-plot_std*8, -plot_std*7, -plot_std*6, -plot_std*5, -plot_std*4, -plot_std*3, -plot_std*2, -plot_std*1, 0, plot_std*1, plot_std*2, plot_std*3]);
xticklabels(num2cell(-8:3));



set(gca,'box','off');
set(gca,'tickdir','out');
hold off