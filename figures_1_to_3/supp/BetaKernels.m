clc;clear all;close all;
addpath('..\lib\')
%%
cPath = 'X:\Widefield'; animals = {'mSM63','mSM64','mSM65','mSM66'};
glmPath = 'X:\Widefield\glm_hmm_models\global_model_map.mat';
mintrialnum = 20; %the minimum number of trials per state to be included in plotting
sessiondates = getGlobalGLMHMMSessions(glmPath); %get sessions with GLM-HMM data
fileprefix = 'final';

%% Retrain models over different states
for i = 1:length(animals)
    for j = 1:length(sessiondates{i})
        fprintf('\nRunning for %s, %s.\n\n',animals{i},sessiondates{i}{j});
        runRidge_overStates(cPath,animals{i},sessiondates{i}{j},glmPath);
    end
end
%% now look at some beta kernels for all sessions

% REG = 'bhvVideo';
REG = 'whiskAnalog';
% REG = 'piezomoveAnalog';
% REG = 'bodyAnalog';
% REG = 'noseAnalog';
% REG = 'faceAnalog';
% REG = 'Move';
% REG = 'slowPupil';

all_eng = {};
all_bias = {};

count = 1;
for i = 1:length(animals)
    for j = 1:length(sessiondates{i})
        fprintf('\nRunning for %s, %s.\n\n',animals{i},sessiondates{i}{j});

        [regLabels, kernels, U] = returnAnalogBetaKernels(animals{i},sessiondates{i}{j},[fileprefix 'fullA.mat'], REG);

        if size(regLabels,2) == 0
            continue %skip, there is no encoding model
        end

        transParams = loadTransParams(cPath, animals{i}, sessiondates{i}{j});
        regind = find(strcmpi(regLabels,REG));
        regbetas = kernels;
        regbetas = unSVDalign2allen(regbetas',U,transParams,[],false);
        %regbetas(isnan(regbetas)) = 0;
        all_eng{count} = regbetas;

        [regLabels, kernels, U] = returnAnalogBetaKernels(animals{i},sessiondates{i}{j},[fileprefix 'fullB.mat'], REG);
        transParams = loadTransParams(cPath, animals{i}, sessiondates{i}{j});
        regind = find(strcmpi(regLabels,REG));

        if isempty(regind)
            fprintf('no regressor in this model')
            continue;
        end

        regbetas = kernels;
        regbetas = unSVDalign2allen(regbetas',U,transParams,[],false);
        %regbetas(isnan(regbetas)) = 0;
        all_bias{count} = regbetas;
        animal_id(count) = i;
        count = count+1;
    end
end

all_eng = cat(4,all_eng{:});
all_bias = cat(4,all_bias{:});

si = 35;
eng_trimmed = all_eng(:,:,:,si);
bias_trimmed = all_bias(:,:,:,si);

%%
clims = [-.0005, .0005];
clims = [-.03 .03];
cmap = 'colormap_blueblackred';

combo = horzcat(eng_trimmed, bias_trimmed, eng_trimmed - bias_trimmed);

writerObj = VideoWriter(['C:\Data\churchland\state_manuscript_new_figs\raw_figures\beta_weights\' REG '_Session' num2str(si) '.avi']);
writerObj.FrameRate = 2;
open(writerObj);

for i = 1:size(combo,3)
    fig = figure;
    title = ['Regressor: ' REG '. Displaying dimension ' num2str(i)];
    plotHeatmap(combo(:,:,i), clims, title, 'Beta weight', cmap, 12);
    writeVideo(writerObj,getframe(gcf));
    close(fig)
end
close(writerObj)

%% pub quality plots for beta weight heatmaps
cl = .0005;
clims = [-.08 .08];ind2plot = 2;
clims = [-.2, .2];ind2plot = 1;
clims = [-.05, .05];ind2plot = 4;
%clims = [-.04 .04];ind2plot = 5;
figure
subplot(1,3,1)
plotHeatmap(eng_trimmed(:,:,ind2plot), clims, ['Engaged... '], 'Beta weight', cmap, 12);
%plotHeatmap(eng_trimmed(:,:), clims, ['Engaged... frame ' num2str(frame2plot)], 'Beta weight', cmap, 12);
subplot(1,3,2)
plotHeatmap(bias_trimmed(:,:,ind2plot), clims, ['Biased... '], 'Beta weight', cmap, 12);
%plotHeatmap(bias_trimmed(:,:), clims, ['Biased... frame ' num2str(frame2plot)], 'Beta weight', cmap, 12);
subplot(1,3,3)
plotHeatmap(eng_trimmed(:,:,ind2plot) - bias_trimmed(:,:,ind2plot), clims, ['Biased... '], 'Beta weight', cmap, 12);

%% get ROI
mask = createCircleMask(size(all_eng,[1,2]), 415, 286 ,30); %

roi_eng = mask .* all_eng(:,:,ind2plot,:);
roi_eng(roi_eng == 0) = NaN;
roi_eng = squeeze(squeeze(mean(roi_eng, [1,2], 'omitnan')));

roi_bias = mask .* all_bias(:,:,ind2plot,:);
roi_bias(roi_bias == 0) = NaN;
roi_bias = squeeze(squeeze(mean(roi_bias, [1,2], 'omitnan')));

figure; hold on;
time = (1:size(roi_eng,1)) / 15; % divide by fs to get time
%plot(time, roi_eng,'LineWidth',2);
%plot(time, roi_bias,'LineWidth',2);
stdshade(roi_eng',.2,'red',time,6,[1],[]);
stdshade(roi_bias',.2,'blue',time,6,[1],[]);

% now plot the mask that was used
figure; hold on;
mask = double(mask);
mask(mask==0) = NaN;
figure
plotHeatmap(eng_trimmed, clims, [], 'Beta weight', cmap, 12);
hold on
plotHeatmap(mask, clims, [], 'Beta weight', cmap, 12);
title('')

%% stats
scores = double([roi_eng; roi_bias]); 
grouped_animals = nominal([animal_id, animal_id]);
is_engaged = nominal([ones(size(roi_eng,1),1)', zeros(size(roi_bias,1),1)']);
pair_id = nominal([1:1:size(roi_eng,1), 1:1:size(roi_bias,1)]);

data = table(grouped_animals', scores, is_engaged', pair_id', 'VariableNames', {'Animals' ,'Score','IsEngaged','PairID'});

formula = 'Score ~ IsEngaged + (IsEngaged|Animals) + (1|PairID)'; % Test if the overall mean is different from 0, accounting for group variability
lme = fitlme(data, formula);
disp(lme)
%lme.Coefficients.pValue


