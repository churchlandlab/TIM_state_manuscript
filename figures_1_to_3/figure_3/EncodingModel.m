clc;clear all;close all;
%% Get the animals and sessions
addpath('..\lib\')


cPath = 'X:\Widefield'; animals = {'mSM63','mSM64','mSM65','mSM66'};

glmPath = 'X:\Widefield\glm_hmm_models\map_all_subjects_targrate.mat';
sessiondates = getGlobalGLMHMMSessions(glmPath); %get sessions with GLM-HMM data

fileprefix = 'final'; %map_all_subjects_targrate



%% Retrain models over different states

%runRidge_overStates(cPath,'CSP22','23-Jun-2020',glmPath);

for i = 1:length(animals)
    parfor j = 1:length(sessiondates{i})
        fprintf('\nRunning for %s, %s.\n\n',animals{i},sessiondates{i}{j});

        runRidge_overStates(cPath,animals{i},sessiondates{i}{j},glmPath, fileprefix);

    end
end


%% get the data
counter = 1;
for i = 1:length(animals)
    for j = 2:length(sessiondates{i})
        sessionname{counter} = [animals{i}, '_', sessiondates{i}{j}];
        fulla(counter,:) = returnVarianceMovie(cPath, animals{i},sessiondates{i}{j}, [fileprefix 'fullA.mat']);
        sponta(counter,:) = returnVarianceMovie(cPath, animals{i},sessiondates{i}{j}, [fileprefix 'spontA.mat']);
        opa(counter,:) = returnVarianceMovie(cPath, animals{i},sessiondates{i}{j}, [fileprefix 'operantA.mat']);
        taskvara(counter,:) = returnVarianceMovie(cPath, animals{i},sessiondates{i}{j}, [fileprefix 'taskA.mat']);
        nosponta(counter,:) = returnVarianceMovie(cPath, animals{i},sessiondates{i}{j}, [fileprefix 'nospontA.mat']);
        noopa(counter,:) = returnVarianceMovie(cPath, animals{i},sessiondates{i}{j}, [fileprefix 'nooperantA.mat']);
        notaskvara(counter,:) = returnVarianceMovie(cPath, animals{i},sessiondates{i}{j}, [fileprefix 'notaskA.mat']);

        fullb(counter,:) = returnVarianceMovie(cPath, animals{i},sessiondates{i}{j}, [fileprefix 'fullB.mat']);
        spontb(counter,:) = returnVarianceMovie(cPath, animals{i},sessiondates{i}{j}, [fileprefix 'spontB.mat']);
        opb(counter,:) = returnVarianceMovie(cPath, animals{i},sessiondates{i}{j}, [fileprefix 'operantB.mat']);
        taskvarb(counter,:) = returnVarianceMovie(cPath, animals{i},sessiondates{i}{j}, [fileprefix 'taskB.mat']);
        nospontb(counter,:) = returnVarianceMovie(cPath, animals{i},sessiondates{i}{j}, [fileprefix 'nospontB.mat']);
        noopb(counter,:) = returnVarianceMovie(cPath, animals{i},sessiondates{i}{j}, [fileprefix 'nooperantB.mat']);
        notaskvarb(counter,:) = returnVarianceMovie(cPath, animals{i},sessiondates{i}{j}, [fileprefix 'notaskB.mat']);

        counter = counter + 1;
        fprintf('\ncounter is %i\n',counter);
    end
end

dsponta = fulla-nosponta;
dspontb = fullb-nospontb;
dtaska = fulla-notaskvara;
dtaskb = fullb-notaskvarb;
dopa = fulla-noopa;
dopb = fullb-noopb;

%% plotting
time = linspace(0,5,size(fulla,2));
time = time-time(30);

figure; hold on;
title('Full Model')
%plot(fulla,'r')
%plot(fullb,'b')
stdshade(fulla,.2,'red',time,6,[30],[]);
stdshade(fullb,.2,'blue',time,6,[30],[]);
ylabel('cvR^2');
xlabel('Time from handle grab (s)')
legend({'','','Engaged','','','Disengaged','','','',''})

figure; hold on;
title('Task Variables')
stdshade(taskvara,.2,'red',time,6,[30],[]);
stdshade(taskvarb,.2,'blue',time,6,[30],[]);
ylabel('cvR^2');
xlabel('Time from handle grab (s)')
legend({'','','Engaged','','','Disengaged','','','',''})

figure; hold on;
title('Instructed Movement Variables')
stdshade(opa,.2,'red',time,6,[30],[]);
stdshade(opb,.2,'blue',time,6,[30],[]);
ylabel('cvR^2');
xlabel('Time from handle grab (s)')
legend({'','','Engaged','','','Disengaged','','','',''})

figure; hold on;
title('Uninstructed Movement Variables')
stdshade(sponta,.2,'red',time,6,[30],[]);
stdshade(spontb,.2,'blue',time,6,[30],[]);
ylabel('cvR^2');
xlabel('Time from handle grab (s)')
legend({'','','Engaged','','','Disengaged','','','',''})



figure; hold on;
title('Task Variables')
stdshade(dtaska,.2,'red',time,6,[30],[]);
stdshade(dtaskb,.2,'blue',time,6,[30],[]);
ylabel('deltaR^2');
xlabel('Time from handle grab (s)')
legend({'','','Engaged','','','Disengaged','','','',''})

figure; hold on;
title('Instructed Movement Variables')
stdshade(dopa,.2,'red',time,6,[30],[]);
stdshade(dopb,.2,'blue',time,6,[30],[]);
ylabel('deltaR^2');
xlabel('Time from handle grab (s)')
legend({'','','Engaged','','','Disengaged','','','',''})

figure; hold on;
title('Uninstructed Movement Variables')
stdshade(dsponta,.2,'red',time,6,[30],[]);
stdshade(dspontb,.2,'blue',time,6,[30],[]);
ylabel('deltaR^2');
xlabel('Time from handle grab (s)')
legend({'','','Engaged','','','Disengaged','','','',''})

%exportgraphics(gcf,'C:\Data\churchland\PowerpointsPostersPresentations\SFN2022/FridayUpdate\encodingmodel\fullcvr.pdf');

%% get the data - but now with better alignment

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

        fulla(counter,:) = returnRealignedVarianceMovie(cPath, animals{i},sessiondates{i}{j}, [fileprefix 'fullA.mat'], segIdx, NFRAMES);
        sponta(counter,:) = returnRealignedVarianceMovie(cPath, animals{i},sessiondates{i}{j}, [fileprefix 'spontA.mat'], segIdx, NFRAMES);
        opa(counter,:) = returnRealignedVarianceMovie(cPath, animals{i},sessiondates{i}{j}, [fileprefix 'operantA.mat'], segIdx, NFRAMES);
        taskvara(counter,:) = returnRealignedVarianceMovie(cPath, animals{i},sessiondates{i}{j}, [fileprefix 'taskA.mat'], segIdx, NFRAMES);
        nosponta(counter,:) = returnRealignedVarianceMovie(cPath, animals{i},sessiondates{i}{j}, [fileprefix 'nospontA.mat'], segIdx, NFRAMES);
        noopa(counter,:) = returnRealignedVarianceMovie(cPath, animals{i},sessiondates{i}{j}, [fileprefix 'nooperantA.mat'], segIdx, NFRAMES);
        notaskvara(counter,:) = returnRealignedVarianceMovie(cPath, animals{i},sessiondates{i}{j}, [fileprefix 'notaskA.mat'], segIdx, NFRAMES);

        fullb(counter,:) = returnRealignedVarianceMovie(cPath, animals{i},sessiondates{i}{j}, [fileprefix 'fullB.mat'], segIdx, NFRAMES);
        spontb(counter,:) = returnRealignedVarianceMovie(cPath, animals{i},sessiondates{i}{j}, [fileprefix 'spontB.mat'], segIdx, NFRAMES);
        opb(counter,:) = returnRealignedVarianceMovie(cPath, animals{i},sessiondates{i}{j}, [fileprefix 'operantB.mat'], segIdx, NFRAMES);
        taskvarb(counter,:) = returnRealignedVarianceMovie(cPath, animals{i},sessiondates{i}{j}, [fileprefix 'taskB.mat'], segIdx, NFRAMES);
        nospontb(counter,:) = returnRealignedVarianceMovie(cPath, animals{i},sessiondates{i}{j}, [fileprefix 'nospontB.mat'], segIdx, NFRAMES);
        noopb(counter,:) = returnRealignedVarianceMovie(cPath, animals{i},sessiondates{i}{j}, [fileprefix 'nooperantB.mat'], segIdx, NFRAMES);
        notaskvarb(counter,:) = returnRealignedVarianceMovie(cPath, animals{i},sessiondates{i}{j}, [fileprefix 'notaskB.mat'], segIdx, NFRAMES);

        counter = counter + 1;
        fprintf('\ncounter is %i\n',counter);
    end
end

dsponta = fulla-nosponta;
dspontb = fullb-nospontb;
dtaska = fulla-notaskvara;
dtaskb = fullb-notaskvarb;
dopa = fulla-noopa;
dopb = fullb-noopb;

%% plotting
fs = 15;
time = 0:1/fs:(size(fulla,2)-1)/fs;
naninds = cumsum(floor(segIdx * fs));
naninds = naninds(1:end-1);


figure; hold on;
title('Full Model')
%plot(fulla,'r')
%plot(fullb,'b')
stdshade(fulla,.2,'red',time,6,naninds,[]);
stdshade(fullb,.2,'blue',time,6,naninds,[]);
ylabel('cvR^2');
xlabel('Time (s)')
legend({'','','','','','Engaged','','','','','','Disengaged'})

figure; hold on;
title('Task Variables')
stdshade(taskvara,.2,'red',time,6,naninds,[]);
stdshade(taskvarb,.2,'blue',time,6,naninds,[]);
ylabel('cvR^2');
xlabel('Time (s)')
legend({'','','','','','Engaged','','','','','','Disengaged'})

figure; hold on;
title('Instructed Movement Variables')
stdshade(opa,.2,'red',time,6,naninds,[]);
stdshade(opb,.2,'blue',time,6,naninds,[]);
ylabel('cvR^2');
xlabel('Time (s)')
legend({'','','','','','Engaged','','','','','','Disengaged'})

figure; hold on;
title('Uninstructed Movement Variables')
stdshade(sponta,.2,'red',time,6,naninds,[]);
stdshade(spontb,.2,'blue',time,6,naninds,[]);
ylabel('cvR^2');
xlabel('Time (s)')
legend({'','','','','','Engaged','','','','','','Disengaged'})



figure; hold on;
title('Task Variables')
stdshade(dtaska,.2,'red',time,6,naninds,[]);
stdshade(dtaskb,.2,'blue',time,6,naninds,[]);
ylabel('deltaR^2');
xlabel('Time (s)')
legend({'','','','','','Engaged','','','','','','Disengaged'})

figure; hold on;
title('Instructed Movement Variables')
stdshade(dopa,.2,'red',time,6,naninds,[]);
stdshade(dopb,.2,'blue',time,6,naninds,[]);
ylabel('deltaR^2');
xlabel('Time (s)')
legend({'','','','','','Engaged','','','','','','Disengaged'})

figure; hold on;
title('Uninstructed Movement Variables')
stdshade(dsponta,.2,'red',time,6,naninds,[]);
stdshade(dspontb,.2,'blue',time,6,naninds,[]);
ylabel('deltaR^2');
xlabel('Time (s)')
legend({'','','','','','Engaged','','','','','','Disengaged'})

%% stats
stim_a = fulla(:,naninds(2):naninds(3)-1);
%stim_a = reshape(stim_a,numel(stim_a),[]);
stim_a = nanmean(stim_a,2);
stim_b = fullb(:,naninds(2):naninds(3)-1);
%stim_b = reshape(stim_b,numel(stim_b),[]);
stim_b = nanmean(stim_b,2);
[h,p] = ttest(stim_a,stim_b)

delay_a = fulla(:,naninds(3):naninds(4)-1);
delay_a = nanmean(delay_a,2);
delay_b = fullb(:,naninds(3):naninds(4)-1);
delay_b = nanmean(delay_b,2);
[h,p] = ttest(delay_a,delay_b)

delay_a = fulla(:,naninds(1):naninds(2));
delay_a = nanmean(delay_a,2);
delay_b = fullb(:,naninds(1):naninds(2));
delay_b = nanmean(delay_b,2);
[h,p] = ttest(delay_a,delay_b)



%% plot a variance map

counter = 1;
for i = 1:length(animals)
    for j = 1:length(sessiondates{i})

        fulla(counter,:,:) = returnVarianceMap(cPath, animals{i},sessiondates{i}{j}, [fileprefix 'fullA.mat']);
        sponta(counter,:,:) = returnVarianceMap(cPath, animals{i},sessiondates{i}{j}, [fileprefix 'spontA.mat']);
        opa(counter,:,:) = returnVarianceMap(cPath, animals{i},sessiondates{i}{j}, [fileprefix 'operantA.mat']);
        taskvara(counter,:,:) = returnVarianceMap(cPath, animals{i},sessiondates{i}{j}, [fileprefix 'taskA.mat']);
        nosponta(counter,:,:) = returnVarianceMap(cPath, animals{i},sessiondates{i}{j}, [fileprefix 'nospontA.mat']);
        noopa(counter,:,:) = returnVarianceMap(cPath, animals{i},sessiondates{i}{j}, [fileprefix 'nooperantA.mat']);
        notaskvara(counter,:,:) = returnVarianceMap(cPath, animals{i},sessiondates{i}{j}, [fileprefix 'notaskA.mat']);

        fullb(counter,:,:) = returnVarianceMap(cPath, animals{i},sessiondates{i}{j}, [fileprefix 'fullB.mat']);
        spontb(counter,:,:) = returnVarianceMap(cPath, animals{i},sessiondates{i}{j}, [fileprefix 'spontB.mat']);
        opb(counter,:,:) = returnVarianceMap(cPath, animals{i},sessiondates{i}{j}, [fileprefix 'operantB.mat']);
        taskvarb(counter,:,:) = returnVarianceMap(cPath, animals{i},sessiondates{i}{j}, [fileprefix 'taskB.mat']);
        nospontb(counter,:,:) = returnVarianceMap(cPath, animals{i},sessiondates{i}{j}, [fileprefix 'nospontB.mat']);
        noopb(counter,:,:) = returnVarianceMap(cPath, animals{i},sessiondates{i}{j}, [fileprefix 'nooperantB.mat']);
        notaskvarb(counter,:,:) = returnVarianceMap(cPath, animals{i},sessiondates{i}{j}, [fileprefix 'notaskB.mat']);

        counter = counter + 1;
        fprintf('\ncounter is %i\n',counter);
    end
end

dsponta = fulla-nosponta;
dspontb = fullb-nospontb;
dtaska = fulla-notaskvara;
dtaskb = fullb-notaskvarb;
dopa = fulla-noopa;
dopb = fullb-noopb;
%%
clims = [0 .8];
clims2 = [-.15 .15];

a = squeeze(nanmean(fulla,1));
b = squeeze(nanmean(fullb,1));

figure;
mapImg = imshow(a, clims);
colormap(mapImg.Parent,'inferno'); axis image;
set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
hcb = colorbar;
%hcb.Title.String = 'cvR^2';

figure;
mapImg = imshow(a, clims);
colormap(mapImg.Parent,'inferno'); axis image;
set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
hcb = colorbar;
%hcb.Title.String = 'cvR^2';

figure;
mapImg = imshow(a-b, clims2);
colormap(mapImg.Parent,'colormap_blueblackred'); axis image;
set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
hcb = colorbar;
%hcb.Title.String = 'cvR^2';



