clc;clear all;close all;
%% Get the animals and sessions
addpath('..\lib\')

cPath = 'X:\Widefield'; animals = {'mSM63','mSM64','mSM65','mSM66'};
glmPath = 'X:\Widefield\glm_hmm_models\global_model_map.mat';
sessiondates = getGlobalGLMHMMSessions(glmPath); %get sessions with GLM-HMM data


fileprefix = 'task_aligned_unaligned';

%% Retrain models over different states

for i = 1:length(animals)
    for j = 1:length(sessiondates{i})
        fprintf('\nRunning for %s, %s.\n\n',animals{i},sessiondates{i}{j});
        %runRidge_taskIndependentVariance(cPath,animals{i},sessiondates{i}{j},glmPath, fileprefix);
        runRidge_taskIndependentVariance_grouped(cPath,animals{i},sessiondates{i}{j},glmPath, fileprefix);
    end
end

%% get the data

allmotorlabels = {{'lGrab','lGrabRel','rGrab','rGrabRel'},{'lLick','rLick'},{'piezoAnalog','piezoDigital','piezoMoveAnalog','piezoMoveDigital','piezoMoveHiDigital'},{'whiskAnalog','whiskDigital','whiskHiDigital'},{'noseAnalog','noseDigital','noseHiDigital'},{'fastPupilAnalog','fastPupilDigital','fastPupilHiDigital','slowPupil'},{'faceAnalog','faceDigital','faceHiDigital'},{'bodyAnalog','bodyDigital','bodyHiDigital'},{'Move'},{'bhvVideo'}};
for i = 1:length(allmotorlabels)
    Afilenames{i} = [fileprefix strjoin(allmotorlabels{i}, '_') 'plustaskA'];
    Bfilenames{i} = [fileprefix strjoin(allmotorlabels{i}, '_') 'plustaskB'];
    singleVarAfilenames{i} = [fileprefix strjoin(allmotorlabels{i}, '_') 'A'];
    singleVarBfilenames{i} = [fileprefix strjoin(allmotorlabels{i}, '_') 'B'];
    figureTitles{i} = strjoin(allmotorlabels{i}, '_');
end

taskIndependentVarianceA = []; taskDependentVarianceA = [];
taskIndependentVarianceB = []; taskDependentVarianceB = [];
for k = 1:length(Afilenames)
    counter = 1;
    for i = 1:length(animals)
        for j = 1:length(sessiondates{i})
            a(counter,:) = returnVariance(cPath, animals{i},sessiondates{i}{j}, Afilenames{k});
            b(counter,:) = returnVariance(cPath, animals{i},sessiondates{i}{j}, Bfilenames{k});
            singleVara(counter,:) = returnVariance(cPath, animals{i},sessiondates{i}{j}, singleVarAfilenames{k});
            singleVarb(counter,:) = returnVariance(cPath, animals{i},sessiondates{i}{j}, singleVarBfilenames{k});
            taska(counter,:) = returnVariance(cPath, animals{i},sessiondates{i}{j}, [fileprefix 'taskA']);
            taskb(counter,:) = returnVariance(cPath, animals{i},sessiondates{i}{j}, [fileprefix 'taskB']);

            counter = counter + 1;
            fprintf('\ncounter is %i\n',counter);
        end
    end
    taskIndependentVarianceA = [taskIndependentVarianceA, a - taska];
    taskIndependentVarianceB = [taskIndependentVarianceB, b - taskb];
    taskDependentVarianceA = [taskDependentVarianceA, singleVara - (a - taska)];
    taskDependentVarianceB = [taskDependentVarianceB, singleVarb - (b - taskb)];
end

%% plotting

figureTitles = {'handles','licks','piezo','whisk','nose','pupil','face','body','videoME','video'};

animalindex = [];
for i = 1:length(animals)
    animalindex = [animalindex; repmat(i, length(sessiondates{i}),1)];
end

onevec = ones(size(taskDependentVarianceA(:,1),1),1);
twovec = onevec.*2;
threevec = onevec.*3;
fourvec = onevec.*4;

for i = 1:size(taskDependentVarianceA,2)
    figure('Position',[500 500 900 400]); hold on;
    title(figureTitles{i})
    scatter(onevec,taskDependentVarianceA(:,i),'b')
    scatter(twovec,taskDependentVarianceB(:,i),'r')
    scatter(threevec,taskIndependentVarianceA(:,i),'b')
    scatter(fourvec,taskIndependentVarianceB(:,i),'r')

    parallelcoords([taskDependentVarianceA(:,i),taskDependentVarianceB(:,i),taskIndependentVarianceA(:,i),taskIndependentVarianceB(:,i)])
    ylabel('R^2')
    xlim([.5 4.5])
    xticks([1.5 3.5])
    xticklabels({'Task Dependent Variance','Task Independent Variance'})
    legend('Engaged','Disengaged','Location','southeast')


end
%% alternative plotting

for i = 1:size(taskDependentVarianceA,2)
    figure('Position',[600 500 900 400]); hold on;
    title(figureTitles{i})
    scatter(onevec,taskDependentVarianceA(:,i),'b')
    scatter(twovec,taskIndependentVarianceA(:,i),'r')
    scatter(threevec,taskDependentVarianceB(:,i),'b')
    scatter(fourvec,taskIndependentVarianceB(:,i),'r')

    parallelcoords([taskDependentVarianceA(:,i),taskIndependentVarianceA(:,i),taskDependentVarianceB(:,i),taskIndependentVarianceB(:,i)])
    ylabel('R^2')
    xlim([.5 4.5])
    xticks([1.5 3.5])
    xticklabels({'Engaged','Disengaged'})
    legend('Task Dependent Variance','Task Independent Variance','Location','southeast')


end

%% Plotting over time

allmotorlabels = {{'lGrab','lGrabRel','rGrab','rGrabRel'},{'lLick','rLick'},{'piezoAnalog','piezoDigital','piezoMoveAnalog','piezoMoveDigital','piezoMoveHiDigital'},{'whiskAnalog','whiskDigital','whiskHiDigital'},{'noseAnalog','noseDigital','noseHiDigital'},{'fastPupilAnalog','fastPupilDigital','fastPupilHiDigital','slowPupil'},{'faceAnalog','faceDigital','faceHiDigital'},{'bodyAnalog','bodyDigital','bodyHiDigital'},{'Move'},{'bhvVideo'}};
for i = 1:length(allmotorlabels)
    Afilenames{i} = [fileprefix strjoin(allmotorlabels{i}, '_') 'plustaskA'];
    Bfilenames{i} = [fileprefix strjoin(allmotorlabels{i}, '_') 'plustaskB'];
    singleVarAfilenames{i} = [fileprefix strjoin(allmotorlabels{i}, '_') 'A'];
    singleVarBfilenames{i} = [fileprefix strjoin(allmotorlabels{i}, '_') 'B'];
    figureTitles{i} = strjoin(allmotorlabels{i}, '_');
end

a = [];b = [];singleVara=[];singleVarb=[];taska=[];taskb=[];
for k = 1:length(Afilenames)
    taskIndependentVarianceA = []; taskDependentVarianceA = [];
    taskIndependentVarianceB = []; taskDependentVarianceB = [];
    counter = 1;
    for i = 1:length(animals)
        for j = 1:length(sessiondates{i})

            if sum(ismember(animals{i},'mSM')) == 3 %mSM Mice
                segIdx = [1 0.5 1.00 0.75 .75]; %[baseline, handle, stim, delay, response] maximal duration of each segment in seconds, use this for EMX mice
            elseif sum(ismember(animals{i},'CSP')) == 3 %CSP Mice
                segIdx = [1 0.5 1.00 0.4 .75]; %testing
            end
            a(counter,:) = returnVarianceMovie(cPath, animals{i},sessiondates{i}{j}, Afilenames{k});
            b(counter,:) = returnVarianceMovie(cPath, animals{i},sessiondates{i}{j}, Bfilenames{k});
            singleVara(counter,:) = returnVarianceMovie(cPath, animals{i},sessiondates{i}{j}, singleVarAfilenames{k});
            singleVarb(counter,:) = returnVarianceMovie(cPath, animals{i},sessiondates{i}{j}, singleVarBfilenames{k});
            taska(counter,:) = returnVarianceMovie(cPath, animals{i},sessiondates{i}{j}, [fileprefix 'taskA']);
            taskb(counter,:) = returnVarianceMovie(cPath, animals{i},sessiondates{i}{j}, [fileprefix 'taskB']);


            counter = counter + 1;
            fprintf('\ncounter is %i\n',counter);
        end
    end
    taskIndependentVarianceA = a - taska;
    taskIndependentVarianceB = b - taskb;
    taskDependentVarianceA = singleVara - (a - taska);
    taskDependentVarianceB = singleVarb - (b - taskb);

    handleframe = 30;
    figure('Position',[600 500 900 400]); hold on;
    title(figureTitles{k})
    stdshade(taskDependentVarianceA,.2,'red',[],6,[handleframe],[])
    stdshade(taskIndependentVarianceA,.2,'blue',[],6,[handleframe],[])
    stdshade(taskDependentVarianceB,.2,'green',[],6,[handleframe],[])
    stdshade(taskIndependentVarianceB,.2,'magenta',[],6,[handleframe],[])

    ylabel('R^2');
    leg = {'','','Task Dependent - engaged','','','Task Independent - engaged','','','Task Dependent - disengaged','','','Task Independent - disengaged'};
    %legend(leg,'Location','southwest');
    legend(leg);


end

%% Plotting over time - better alignment

NFRAMES = 75;

allmotorlabels = {{'lGrab','lGrabRel','rGrab','rGrabRel'},{'lLick','rLick'},{'piezoAnalog','piezoDigital','piezoMoveAnalog','piezoMoveDigital','piezoMoveHiDigital'},{'whiskAnalog','whiskDigital','whiskHiDigital'},{'noseAnalog','noseDigital','noseHiDigital'},{'fastPupilAnalog','fastPupilDigital','fastPupilHiDigital','slowPupil'},{'faceAnalog','faceDigital','faceHiDigital'},{'bodyAnalog','bodyDigital','bodyHiDigital'},{'Move'},{'bhvVideo'}};
for i = 1:length(allmotorlabels)
    Afilenames{i} = [fileprefix strjoin(allmotorlabels{i}, '_') 'plustaskA'];
    Bfilenames{i} = [fileprefix strjoin(allmotorlabels{i}, '_') 'plustaskB'];
    singleVarAfilenames{i} = [fileprefix strjoin(allmotorlabels{i}, '_') 'A'];
    singleVarBfilenames{i} = [fileprefix strjoin(allmotorlabels{i}, '_') 'B'];
    figureTitles{i} = strjoin(allmotorlabels{i}, '_');
end

a = [];b = [];singleVara=[];singleVarb=[];taska=[];taskb=[];
for k = 1:length(Afilenames)
    taskIndependentVarianceA = []; taskDependentVarianceA = [];
    taskIndependentVarianceB = []; taskDependentVarianceB = [];
    counter = 1;
    for i = 1:length(animals)
        for j = 1:length(sessiondates{i})

            if sum(ismember(animals{i},'mSM')) == 3 %mSM Mice
                segIdx = [1 0.5 1.00 0.75 .75]; %[baseline, handle, stim, delay, response] maximal duration of each segment in seconds, use this for EMX mice
            elseif sum(ismember(animals{i},'CSP')) == 3 %CSP Mice
                segIdx = [1 0.5 1.00 0.4 .75]; %testing
            end

            a(counter,:) = returnRealignedVarianceMovie(cPath, animals{i},sessiondates{i}{j}, Afilenames{k}, segIdx, NFRAMES);
            b(counter,:) = returnRealignedVarianceMovie(cPath, animals{i},sessiondates{i}{j}, Bfilenames{k}, segIdx, NFRAMES);
            singleVara(counter,:) = returnRealignedVarianceMovie(cPath, animals{i},sessiondates{i}{j}, singleVarAfilenames{k}, segIdx, NFRAMES);
            singleVarb(counter,:) = returnRealignedVarianceMovie(cPath, animals{i},sessiondates{i}{j}, singleVarBfilenames{k}, segIdx, NFRAMES);
            taska(counter,:) = returnRealignedVarianceMovie(cPath, animals{i},sessiondates{i}{j}, [fileprefix 'taskA'], segIdx, NFRAMES);
            taskb(counter,:) = returnRealignedVarianceMovie(cPath, animals{i},sessiondates{i}{j}, [fileprefix 'taskB'], segIdx, NFRAMES);

            counter = counter + 1;
            fprintf('\ncounter is %i\n',counter);
        end
    end
    taskIndependentVarianceA = a - taska;
    taskIndependentVarianceB = b - taskb;
    taskDependentVarianceA = singleVara - (a - taska);
    taskDependentVarianceB = singleVarb - (b - taskb);

    time = 1:size(taskDependentVarianceA,2);
    naninds = cumsum(floor(segIdx * 15));
    naninds = naninds(1:end-1);

    figure('Position',[600 500 900 400]); hold on;
    title(figureTitles{k})
    stdshade(taskDependentVarianceA,.2,'red',time,6,naninds,[])
    stdshade(taskIndependentVarianceA,.2,'blue',time,6,naninds,[])
    stdshade(taskDependentVarianceB,.2,'green',time,6,naninds,[])
    stdshade(taskIndependentVarianceB,.2,'magenta',time,6,naninds,[])

    ylabel('R^2');
    leg = cell(24,1); leg(:) = {''};
    leg{6} = 'Task Dependent - engaged';
    leg{12} = 'Task Independent - engaged';
    leg{18} = 'Task Dependent - disengaged';
    leg{24} = 'Task Independent - disengaged';

    %legend(leg,'Location','southwest');
    legend(leg);


end

