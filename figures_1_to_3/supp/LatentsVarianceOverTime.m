clc;clear all;close all;
addpath('..\lib\')
%% Get the animals and sessions
cPath = 'X:\Widefield'; animals = {'mSM63','mSM64','mSM65','mSM66'};
glmFile = 'X:\Widefield\glm_hmm_models\map_all_subjects_targrate.mat';
sessiondates = getGlobalGLMHMMSessions(glmFile); %get sessions with GLM-HMM data
sessiondates{3} = sessiondates{3}(1:11); % drop mSM65 16-Jul-2018 because it could not be aligned to allen atlas
method = 'cutoff';
mintrialnum = 20;
dualcase = 'reward';
%dualcase = 'choice_equal';
fsize = 29;
clims = {[-.0003 .0003],[-.0003 .0003]}; %for variance

%% get variance over time of the latents (SVD temporal components)
A_latents = []; B_latents = []; animal_inds = [];
count=1;
for j = 1:length(animals)
    for k = 1:length(sessiondates{j})
        fprintf('\n%s, %s',animals{j},sessiondates{j}{k})
        [aa,inds{1},inds{2},~,pstate,~] = getStateInds(cPath,animals{j},sessiondates{j}{k},method,glmFile,dualcase, false);
        p_eng = pstate(1,:);
        p_eng = p_eng(aa);
        if length(inds{1}) < mintrialnum
            fprintf('\nSkipping!\n');
            continue
        end
        [all_latent_vars,eventframes,time] = getLatentsVarianceOverTime(cPath,animals{j},sessiondates{j}{k},inds);
        A_latents(count,:,:) = all_latent_vars{1};
        B_latents(count,:,:) = all_latent_vars{2};
        animal_inds = [animal_inds, j];
        count = count+1;
    end
end
%% plot all sessions
figure; hold on;
dim_to_plot =  1
% Add scatter points
% Calculate means and standard errors
A = squeeze(A_latents(:,dim_to_plot));
B = squeeze(B_latents(:,dim_to_plot));

meanA = mean(A);
meanB = mean(B);
semA = std(A) / sqrt(length(A));
semB = std(B) / sqrt(length(B));

%bp = boxplot([nanmean(tnostate), nanmean(tfull), nanmean(tfullA), nanmean(tfullB)], [1,2,3,4], 'Labels', {'NoState', 'Full', 'Engaged','Disengaged'});
% add means
plot([-.2 .2] + 1, [meanA meanA], 'LineWidth', 2,'Color','black');
plot([-.2 .2] + 2, [meanB meanB], 'LineWidth', 2,'Color','black');
xticks([1,2]);
xticklabels({'Engaged','Disengaged'});

% Add error bars
errorbar(1, meanA, semA, 'k', 'LineWidth', 1.5);
errorbar(2, meanB, semB, 'k', 'LineWidth', 1.5);


sigma = .06;
%cols = {"red","blue","green","magenta"};
cols = ["#D95319", "#EDB120",  "#7E2F8E", "#77AC30"];
plotcols = {};
for i=1:length(cols)
    r = hex2dec(cols{i}(2:3)) / 255;
    g = hex2dec(cols{i}(4:5)) / 255;
    b = hex2dec(cols{i}(6:7)) / 255;
    % Return the RGB triplet
    rgb = [r, g, b];

    valid_inds = animal_inds == i;
    Asmall = A(valid_inds);
    Bsmall = B(valid_inds);
    sc1 = sigma .* randn(sum(valid_inds),1) + 1;
    scatter(sc1, Asmall, [], rgb, 'filled', 'MarkerFaceAlpha', .8);
    sc2 = sigma .* randn(sum(valid_inds),1) + 2;
    scatter(sc2, Bsmall, [], rgb, 'filled', 'MarkerFaceAlpha', .8);
    for j=1:length(Asmall)
        plot([sc1(j), sc2(j)],[Asmall(j), Bsmall(j)],'Color', [0 0 0 0.2]);
    end
end

xlabel(['SVD Dimension', num2str(dim_to_plot)]);
xlim([0,4])
ylabel('Variance over time')
%% Stats 
scores = double([A;B]); 
grouped_animals = nominal([animal_inds, animal_inds]);
is_engaged = nominal([ones(size(A,1),1)', zeros(size(B,1),1)']);
pair_id = nominal([1:1:size(A,1), 1:1:size(B,1)]);

data = table(grouped_animals', scores, is_engaged', pair_id', 'VariableNames', {'Animals' ,'Score','IsEngaged','PairID'});

formula = 'Score ~ IsEngaged + (IsEngaged|Animals) + (1|PairID)'; % Test if the overall mean is different from 0, accounting for group variability
lme = fitlme(data, formula);
disp(lme)


%%
figure; hold on;
for i=1:200
    diff(i,:) = A_latents(:,i) - B_latents(:,i);
end
%stdshade(diff', .2, 'cyan', 1:200, [],[200])
stdshade(A_latents, .2, 'red', 1:200, [],[200])
stdshade(B_latents, .2, 'blue', 1:200, [],[200])
%stdshade(A',.2,'red',time,6,eventframes,[]);

%% do one example session, looking at variance across time
[~,inds{1},inds{2}] = getStateInds(cPath,animals{1},sessiondates{1}{1},method,glmFile,dualcase, false);
[all_latent_vars,eventframes,time] = getLatentsVarianceOverTime(cPath,animals{1},sessiondates{1}{1},inds);
%%
bins = 0:.4:15;
for d=1:5
    f = figure; hold on;
    engaged_trials = all_latent_vars{1}(d,:);
    disengaged_trials = all_latent_vars{2}(d,:);

    bins = 1:.05*max(engaged_trials):max(engaged_trials);
    histogram(engaged_trials,bins,'FaceColor','red')
    histogram(disengaged_trials,bins,'FaceColor','blue')
    xlabel(['Variance in trial: Dimension ', num2str(d)]);
    ylabel('Trials');
    %exportgraphics(f, ['C:\Data\churchland\state_manuscript_new_figs\raw_figures\latents_variance_over_time\individual_trials\dim_' num2str(d) '.pdf'])
end

%% do multiple sessions
A_latents = {}; B_latents = {}; animal_inds = [];
count=1;
for j = 1:length(animals)
    for k = 1:length(sessiondates{j})
        fprintf('\n%s, %s',animals{j},sessiondates{j}{k})
        [aa,inds{1},inds{2},~,pstate,~] = getStateInds(cPath,animals{j},sessiondates{j}{k},method,glmFile,dualcase, false);
        if length(inds{1}) < mintrialnum
            fprintf('\nSkipping!\n');
            continue
        end

        [all_latent_vars,eventframes,time] = getLatentsVarianceOverTime(cPath,animals{j},sessiondates{j}{k},inds);
      
        A_latents{count} = all_latent_vars{1};
        B_latents{count} = all_latent_vars{2};
        animal_inds = [animal_inds, j];
        count = count+1;
    end
end
%%
bins = 0:.4:15;
engaged_trials = cell2mat(cellfun(@(x) mean(x, 2), A_latents, 'UniformOutput', false));
disengaged_trials = cell2mat(cellfun(@(x) mean(x, 2), B_latents, 'UniformOutput', false));

for d=1:5
    figure; hold on;
    histogram(engaged_trials(d,:),bins,'FaceColor','red')
    histogram(disengaged_trials(d,:),bins, 'FaceColor','blue')
    xlabel(['Variance in trial: Dimension ', num2str(d)]);
end
%%
f = figure; hold on;
engaged_trials = cell2mat(cellfun(@(x) mean(x, 2), A_latents, 'UniformOutput', false));
disengaged_trials = cell2mat(cellfun(@(x) mean(x, 2), B_latents, 'UniformOutput', false));

dim_to_plot =  1
% Add scatter points
% Calculate means and standard errors
A = squeeze(engaged_trials(dim_to_plot,:));
B = squeeze(disengaged_trials(dim_to_plot,:));

meanA = mean(A);
meanB = mean(B);
semA = std(A) / sqrt(length(A));
semB = std(B) / sqrt(length(B));

% add means
plot([-.2 .2] + 1, [meanA meanA], 'LineWidth', 2,'Color','black');
plot([-.2 .2] + 2, [meanB meanB], 'LineWidth', 2,'Color','black');
xticks([1,2]);
xticklabels({'Engaged','Disengaged'});

% Add error bars
errorbar(1, meanA, semA, 'k', 'LineWidth', 1.5);
errorbar(2, meanB, semB, 'k', 'LineWidth', 1.5);

sigma = .06;
%cols = {"red","blue","green","magenta"};
cols = ["#D95319", "#EDB120",  "#7E2F8E", "#77AC30"];
plotcols = {};
for i=1:length(cols)
    r = hex2dec(cols{i}(2:3)) / 255;
    g = hex2dec(cols{i}(4:5)) / 255;
    b = hex2dec(cols{i}(6:7)) / 255;
    % Return the RGB triplet
    rgb = [r, g, b];

    valid_inds = animal_inds == i;
    Asmall = A(valid_inds);
    Bsmall = B(valid_inds);
    sc1 = sigma .* randn(sum(valid_inds),1) + 1;
    scatter(sc1, Asmall, [], rgb, 'filled', 'MarkerFaceAlpha', .8);
    sc2 = sigma .* randn(sum(valid_inds),1) + 2;
    scatter(sc2, Bsmall, [], rgb, 'filled', 'MarkerFaceAlpha', .8);
    for j=1:length(Asmall)
        plot([sc1(j), sc2(j)],[Asmall(j), Bsmall(j)],'Color', [0 0 0 0.2]);
    end
end

xlabel(['SVD Dimension', num2str(dim_to_plot)]);
xlim([0,4])
ylabel('Variance over time')
