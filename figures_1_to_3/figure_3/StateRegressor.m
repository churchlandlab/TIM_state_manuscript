clc;clear all;close all;
%% Get the animals and sessions
addpath('..\lib\')

%cPath = 'X:\Widefield'; animals = {'mSM63','mSM64','mSM65','mSM66'};
cPath = 'X:\Widefield'; animals = {'Fez71','Fez72','Fez73','Fez74','Fez75'};

glmPath = 'X:\Widefield\glm_hmm_models\map_all_subjects_targrate.mat';
%glmPath = 'X:\Widefield\glm_hmm_models\fez_4_state.mat';
sessiondates = getGlobalGLMHMMSessions(glmPath); %get sessions with GLM-HMM data
fileprefix = 'state_reg';
%% train the models
for i = 1:length(animals)
    parfor j = 1:length(sessiondates{i})
        fprintf('\nRunning for %s, %s.\n\n',animals{i},sessiondates{i}{j});
        runRidge_stateRegressor(cPath,animals{i},sessiondates{i}{j},glmPath, fileprefix);
        %runRidge_stateRegressor_equalside(cPath,animals{i},sessiondates{i}{j},glmPath, fileprefix);
    end
end
%% get the data
counter = 1;
for i = 1:length(animals)
    for j = 2:length(sessiondates{i})
        sessionname{counter} = [animals{i}, '_', sessiondates{i}{j}];

        full(counter,:) = returnVarianceMovie(cPath, animals{i},sessiondates{i}{j}, [fileprefix 'full.mat']);
        singlestate(counter,:) = returnVarianceMovie(cPath, animals{i},sessiondates{i}{j}, [fileprefix 'singlevar_state.mat']);
        nostate(counter,:) = returnVarianceMovie(cPath, animals{i},sessiondates{i}{j}, [fileprefix 'no_state.mat']);

        counter = counter + 1;
        fprintf('\ncounter is %i\n',counter);
    end
end
dstate = full - nostate;

%% plotting
time = linspace(0,5,size(full,2));
time = time-time(30);

figure; hold on;
title('State regressor')
stdshade(full,.2,'green',time,6,[30],[]);
stdshade(nostate,.2,'cyan',time,6,[30],[]);
ylabel('cvR^2');
xlabel('Time from handle grab (s)')
legend({'','','Full model','','','State regressor removed','','','',''})

figure; hold on;
title('State regressor - single variable')
stdshade(singlestate,.2,'green',time,6,[30],[]);
ylabel('cvR^2');
xlabel('Time from handle grab (s)')
legend({'','','Full model','','','State regressor removed','','','',''})

%% Now with better alignment
NFRAMES = 75;
counter = 1;
for i = 1:length(animals)
    for j = 1:length(sessiondates{i})

        if sum(ismember(animals{i},'mSM')) == 3 %mSM Mice
            segIdx = [1 0.5 1.00 0.75 .75] %[baseline, handle, stim, delay, response] maximal duration of each segment in seconds, use this for EMX mice
        elseif sum(ismember(animals{i},'CSP')) == 3 %CSP Mice
            segIdx = [1 0.5 1.00 0.4 .75] %testing
        elseif sum(ismember(animals{i},'Fez')) == 3
            segIdx = [1 0.5 1.00 0.75 .75]
        end

        full(counter,:) = returnRealignedVarianceMovie(cPath, animals{i},sessiondates{i}{j}, [fileprefix 'full.mat'], segIdx, NFRAMES);
        singlestate(counter,:) = returnRealignedVarianceMovie(cPath, animals{i},sessiondates{i}{j}, [fileprefix 'singlevar_state.mat'], segIdx, NFRAMES);
        nostate(counter,:) = returnRealignedVarianceMovie(cPath, animals{i},sessiondates{i}{j}, [fileprefix 'no_state.mat'], segIdx, NFRAMES);

        counter = counter + 1;
        fprintf('\ncounter is %i\n',counter);
    end
end

dstate = full - nostate;

%% plotting

fs = 15;
time = 0:1/fs:(size(full,2)-1)/fs;
naninds = cumsum(floor(segIdx * fs));
naninds = naninds(1:end-1);

figure; hold on;
title('Full Model - state regressor')
stdshade(full,.2,'green',time,6,naninds,[]);
stdshade(nostate,.2,'cyan',time,6,naninds,[]);
ylabel('cvR^2');
xlabel('Time (s)')
legend({'','','','','','Full model','','','','','','State regressor shuffled'})

figure; hold on;
title('Single variable model - state')
stdshade(singlestate,.2,'green',time,6,naninds,[]);
ylabel('cvR^2');
xlabel('Time (s)')
legend({'','','','','','State regressor','','','','','',''})

%% bar plot
counter = 1;
animal_group = [];
for i = 1:length(animals)
    for j = 1:length(sessiondates{i})
        sessionname{counter} = [animals{i}, '_', sessiondates{i}{j}];
        animal_ind(counter) = i;
        full(counter,:) = returnVarianceMovie(cPath, animals{i},sessiondates{i}{j}, [fileprefix 'full.mat']);
        singlestate(counter,:) = returnVarianceMovie(cPath, animals{i},sessiondates{i}{j}, [fileprefix 'singlevar_state.mat']);
        nostate(counter,:) = returnVarianceMovie(cPath, animals{i},sessiondates{i}{j}, [fileprefix 'no_state.mat']);
        fullA(counter,:) = returnVarianceMovie(cPath, animals{i},sessiondates{i}{j}, [fileprefix 'fullA.mat']);
        fullB(counter,:) = returnVarianceMovie(cPath, animals{i},sessiondates{i}{j}, [fileprefix 'fullB.mat']);

        counter = counter + 1;
        fprintf('\ncounter is %i\n',counter);
    end
end

dstate = full - nostate;

% Combine data into a single matrix and create a grouping variable
tfull = nanmean(full, 2);
tnostate = nanmean(nostate, 2);
tfullA = nanmean(fullA, 2);
tfullB = nanmean(fullB, 2);

% drop nans
animal_ind_1 = animal_ind(~isnan(tfull));
animal_ind_2 = animal_ind(~isnan(tfullA));

tfull = tfull(~isnan(tfull));
tnostate = tnostate(~isnan(tnostate));
tfullA = tfullA(~isnan(tfullA));
tfullB = tfullB(~isnan(tfullB));

figure;
hold on;
%ylim([.2, .4])

set(findobj(gca, 'type', 'line'), 'linew', 2); % Adjust the line width

% Add scatter points
% Calculate means and standard errors
mean_full = mean(tfull);
mean_nostate = mean(tnostate);
meanA = mean(tfullA);
meanB = mean(tfullB);

sem_full = std(tfull) / sqrt(length(tfull));
sem_nostate = std(tnostate) / sqrt(length(tnostate));
sem_A = std(tfullA) / sqrt(length(tfullA));
sem_B = std(tfullB) / sqrt(length(tfullB));

% add means
plot([-.2 .2] + 1, [mean_full mean_full], 'LineWidth', 2,'Color','black');
plot([-.2 .2] + 2, [mean_nostate mean_nostate], 'LineWidth', 2,'Color','black');
plot([-.2 .2] + 3, [meanA meanA], 'LineWidth', 2,'Color','black');
plot([-.2 .2] + 4, [meanB meanB], 'LineWidth', 2,'Color','black');
xticks([1,2,3,4]);
xticklabels({'NoState', 'Full', 'Engaged','Disengaged'});

% Add error bars
errorbar(1, mean_full, sem_full, 'k', 'LineWidth', 1.5);
errorbar(2, mean_nostate, sem_nostate, 'k', 'LineWidth', 1.5);
errorbar(3, meanA, sem_A, 'k', 'LineWidth', 1.5);
errorbar(4, meanB, sem_B, 'k', 'LineWidth', 1.5);

sigma = .06;
cols = ["#D95319", "#EDB120",  "#7E2F8E", "#77AC30",'#FF1493'];
for i = 1:5
    r = hex2dec(cols{i}(2:3)) / 255;
    g = hex2dec(cols{i}(4:5)) / 255;
    b = hex2dec(cols{i}(6:7)) / 255;
    % Return the RGB triplet
    rgb = [r, g, b];

    valid_inds = animal_ind_1 == i;
    sc = sigma .* randn(length(tfull(valid_inds)),1) + 1;
    scatter(sc, tfull(valid_inds), [], rgb, 'filled', 'MarkerFaceAlpha', .8);
    sc = sigma .* randn(length(tnostate(valid_inds)),1) + 2;
    scatter(sc, tnostate(valid_inds), [], rgb, 'filled', 'MarkerFaceAlpha', .8);
    
    valid_inds = animal_ind_2 == i;
    sc1 = sigma .* randn(length(tfullA(valid_inds)),1) + 3;
    scatter(sc1, tfullA(valid_inds), [], rgb, 'filled', 'MarkerFaceAlpha', .8);
    sc2 = sigma .* randn(length(tfullB(valid_inds)),1) + 4;
    scatter(sc2, tfullB(valid_inds), [], rgb, 'filled', 'MarkerFaceAlpha', .8);
    %plot([sc1,sc2]',[tfullA(valid_inds), tfullB(valid_inds)]','Color','black')
    %sc = sigma .* randn(length(tfullB(valid_inds)),1) + 5;
    %scatter(sc, tfullA(valid_inds) - tfullB(valid_inds), [], rgb, 'filled', 'MarkerFaceAlpha', .8);

end

ylabel('cvR^2');
hold off;

%% stats - mixed effects version

scores = [tfull', tnostate'];
grouped_animals = [animal_ind_1, animal_ind_1];
has_state = [ones(length(tfull),1)', zeros(length(tnostate),1)'];
data = table(grouped_animals', has_state', scores', 'VariableNames', {'Animals', 'HasStateRegressor' ,'Score'});
%disp(data);

formula = 'Score ~ 1 + HasStateRegressor + (1|Animals)'; % Test if the overall mean is different from 0, accounting for group variability
lme = fitlme(data, formula);
disp(lme);

%% Now do this for the paired sessions separated by state

scores = [tfullA', tfullB'];
grouped_animals = nominal([animal_ind_2, animal_ind_2]);
is_engaged = nominzal([ones(length(tfullA),1)', zeros(length(tfullB),1)']);
pair_id = nominal([1:1:length(tfullA), 1:1:length(tfullB)]);

data = table(grouped_animals', scores', is_engaged', pair_id', 'VariableNames', {'Animals' ,'Score','IsEngaged','PairID'});
disp(data);

formula = 'Score ~ IsEngaged + (1|Animals) + (1|PairID)'; % Test if the overall mean is different from 0, accounting for group variability
lme = fitlme(data, formula);
disp(lme);