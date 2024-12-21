clc;clear all;close all;
addpath('../lib');
%% Get the animals and sessions

cPath = 'X:\Widefield'; animals = {'mSM63','mSM64','mSM65','mSM66'};
glmFile = 'X:\Widefield\glm_hmm_models\map_all_subjects_targrate.mat';

sessiondates = getGlobalGLMHMMSessions(glmFile); %get sessions with GLM-HMM data
sessiondates{3} = sessiondates{3}(1:11); % drop mSM65 16-Jul-2018 because it could not be aligned to allen atlas

method = 'cutoff';
%method = 'max';
leftChoiceOnly = 1;

mintrialnum = 20; %the minimum number of trials per state to be included in plotting

dualcase = 'reward';
%dualcase = 'choice_equal';
fsize = 29;
clims = {[-.0003 .0003],[-.0003 .0003]}; %for variance
%% Plot avg activity map - individual sessions
inds = {NaN,NaN};
for i = 1:length(animals) %try a few different sessions
    for j = 1:length(sessiondates{i})
        Rec = sessiondates{i}{j};
        fprintf('\nrunning %s on %s\n',animals{i},Rec);
        [~, a, b, ~, ~, SessionData] = getStateInds(cPath,animals{i},Rec,method,glmFile, dualcase, leftChoiceOnly);
        nt = num2str(length(a));
        if length(a) < mintrialnum %skip if too few trials
            out{i,j,:,:} = [];
        else
            out{i,j,:,:} = plotVarianceMap(cPath,animals{i},Rec,{a,b},[animals{i} ' ' Rec ': ' nt ' trials per state'],{'Attentive trials','Bias trials'},clims,fsize,true);
        end
    end
end
clear attend bias
loc = 1; %use this to squish animals and sessions into one dimension
for i = 1:length(animals) %iterate thru animals
    for j = 1:length(sessiondates{i})
        if ~isempty(out{i,j})
            for k = 1:5 %iterate thru trial periods

                attend(loc,k,:,:) = out{i,j}{1,k}; %[animals, trial periods, x, y]
                bias(loc,k,:,:) = out{i,j}{2,k};
            end
            sessionname{loc} = [animals{i}, '_', sessiondates{i}{j}];
            loc = loc + 1;
        end

    end
end
%% 
attendmean = squeeze(mean(attend,1,'omitnan')); %average over animals/sessions
biasmean = squeeze(mean(bias,1,'omitnan')); %average over animals/sessions
combo = cat(4,attendmean,biasmean);

%% plotting the average variance over trials
pltlegend = {'Engaged trials','Bias trials'};
clims = {[0 .00012],[-.00012/2 .00012/2]}; %for variance
fsize = 29;

set(gca,'FontSize',fsize)
%plttitle = 'Activity map averaged over sessions';
plttitle = '';
figure('units','normalized','outerposition',[0 0 1 1],'PaperSize',[40 40])
%figure
for i = 1:2
    subplot(3,5,1+(i-1)*5);
    mapImg = imshow(squeeze(combo(1,:,:,i)), clims{1});
    colormap(mapImg.Parent,'inferno'); axis image;
    if i == 1
        title('Baseline','FontSize',fsize);
    end
    set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
    %hcb = colorbar;
    %hcb.Title.String = 'dF/F';
    ylabel(pltlegend{i},'FontSize',fsize);

    subplot(3,5,2+(i-1)*5);
    mapImg = imshow(squeeze(combo(2,:,:,i)), clims{1});
    colormap(mapImg.Parent,'inferno'); axis image;
    if i == 1
        title('Trial Initiation','FontSize',fsize);
    end
    set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
    %hcb = colorbar;
    %hcb.Title.String = 'dF/F';

    subplot(3,5,3+(i-1)*5);
    mapImg = imshow(squeeze(combo(3,:,:,i)), clims{1});
    colormap(mapImg.Parent,'inferno'); axis image;
    if i == 1
        title('Stimulus','FontSize',fsize);
    end
    set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
    %hcb = colorbar;
    %hcb.Title.String = 'dF/F';

    subplot(3,5,4+(i-1)*5);
    mapImg = imshow(squeeze(combo(4,:,:,i)), clims{1});
    colormap(mapImg.Parent,'inferno'); axis image;
    if i == 1
        title('Delay','FontSize',fsize);
    end
    set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
    %hcb = colorbar;
    %hcb.Title.String = 'dF/F';

    subplot(3,5,5+(i-1)*5);
    mapImg = imshow(squeeze(combo(5,:,:,i)), clims{1});
    colormap(mapImg.Parent,'inferno'); axis image;
    if i == 1
        title('Response','FontSize',fsize);
    end
    set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
    hcb = colorbar;
    hcb.Title.String = 'Variance';
    hcb.Position = hcb.Position + [0.02 0 0 0];
    hcb.FontSize = fsize;
end
mycmap = load('CustomColormap2.mat');
mycmap = mycmap.CustomColormap2;

subplot(3,5,11);
mapImg = imshow(squeeze(combo(1,:,:,1) - combo(1,:,:,2)), clims{2});
colormap(mapImg.Parent,mycmap); axis image; %title('Baseline');
set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
%hcb = colorbar;
%hcb.Title.String = 'dF/F';
ylabel('Difference','FontSize',fsize);

subplot(3,5,12);
mapImg = imshow(squeeze(combo(2,:,:,1) - combo(2,:,:,2)), clims{2});
colormap(mapImg.Parent,mycmap); axis image; %title('Trial Initiation');
set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
%hcb = colorbar;
%hcb.Title.String = 'dF/F';

subplot(3,5,13);
mapImg = imshow(squeeze(combo(3,:,:,1) - combo(3,:,:,2)), clims{2});
colormap(mapImg.Parent,mycmap); axis image; %title('Stimulus');
set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
%hcb = colorbar;
%hcb.Title.String = 'dF/F';

subplot(3,5,14);
mapImg = imshow(squeeze(combo(4,:,:,1) - combo(4,:,:,2)), clims{2});
colormap(mapImg.Parent,mycmap); axis image; %title('Delay');
set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
%hcb = colorbar;
%hcb.Title.String = 'dF/F';

subplot(3,5,15);
mapImg = imshow(squeeze(combo(5,:,:,1) - combo(5,:,:,2)), clims{2});
colormap(mapImg.Parent,mycmap); axis image; %title('Response');
set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
hcb = colorbar;
hcb.Title.String = '\DeltaVariance';
hcb.Position = hcb.Position + [0.01 0 0 0];
hcb.FontSize = fsize;
sgtitle(plttitle);

%% load the allen map, select regions to plot
clear z zname
t = load('C:\Data\churchland\ridgeModel\allenDorsalMapSM.mat');
map = t.dorsalMaps.areaMap;
figure
imagesc(map)
[x,y] = getpts;
x=int64(x);y=int64(y);

for i = 1:length(x)
    z(i) = map(y(i),x(i));
    zname{i} = t.dorsalMaps.labelsSplit(z(i));
end
zname = arrayfun(@string, zname);
fprintf('\nRegion to extract: %s',zname);
close

%% plot variance over time
for i = 1:length(z)
    A = []; B = []; sessionname = {};
    for j = 1:length(animals)
        for k = 1:length(sessiondates{j})
            fprintf('\n%s, %s',animals{j},sessiondates{j}{k})
            [~,inds{1},inds{2}] = getStateInds(cPath,animals{j},sessiondates{j}{k},method,glmFile,dualcase, leftChoiceOnly);

            if length(inds{1}) < mintrialnum
                fprintf('\nSkipping!\n');
                continue
            end

            [temp,eventframes] = plotRegionVar(cPath,animals{j},sessiondates{j}{k},inds,z(i),zname(i),{'','','','','','','Engaged trials','','','','','','','Biased trials'},false,false);
            %exportgraphics(gcf,strjoin(['C:\Data\churchland\PowerpointsPostersPresentations\SFN2022\EMX_individual_psth\' animals{j} sessiondates{j}{k}  zname(i) '.pdf']));
            close gcf;
            A = [A,temp{1}']; B = [B,temp{2}']; %[nframes, ntrials]
            sessionname{end+1} = [animals{j}, '_', sessiondates{j}{k}];
        end
    end
    figure;hold on; title(zname(i))
    time = (0:1:size(A,1)-1) ./ 30; %IMPORTANT, FS IS 15 FOR CSTR MICE
    stdshade(A',.2,'red',time,6,eventframes,[]);
    stdshade(B',.2,'blue',time,6,eventframes,[]);
    legend('','','','','','Engaged trials','','','','','','Biased trials');
    xlabel('Time (s)')
    ylabel('Variance')
    set(gca,'TickDir','out')
    
end




