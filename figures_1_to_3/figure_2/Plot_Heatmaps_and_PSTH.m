% We know that there is a heatmap difference that is state dependent. So lets
% plot this heatmat, then get some PSTH's, and finally, revisit the results
% of the encoding model.
clc;clear all;close all;
addpath('..\lib');
%% Get the animals and sessions


cPath = 'X:\Widefield'; animals = {'mSM63','mSM64','mSM65','mSM66'};

glmFile = 'X:\Widefield\glm_hmm_models\map_all_subjects_targrate.mat';
sessiondates = getGlobalGLMHMMSessions(glmFile); %get sessions with GLM-HMM data
sessiondates{3} = sessiondates{3}(1:11); % drop mSM65 16-Jul-2018 because it could not be aligned to allen atlas

method = 'cutoff';
%method = 'max';
mintrialnum = 20; %the minimum number of trials per state to be included in plotting
%mintrialnum = 25;
dualcase = 'reward';
%dualcase = 'choice_equal';
onlyLeftChoice=0;

fsize = 29;

%clims = {[-.01 .01],[-.005 .005]};
clims = {[-.01 .01],[-.01 .01]}; %for df/f
%% Plot avg activity map - individual sessions
inds = {NaN,NaN};
a_rewarded = 0;b_rewarded = 0;
a_unrewarded = 0; b_unrewarded = 0;
a_right = 0; b_right = 0;
a_left = 0; b_left = 0;
for i = 1:length(animals) %try a few different sessions
    for j = 1:length(sessiondates{i})
        Rec = sessiondates{i}{j};
        fprintf('\nrunning %s on %s\n',animals{i},Rec);
        [~, a, b, ~, ~, SessionData] = getStateInds(cPath,animals{i},Rec,method,glmFile, dualcase, onlyLeftChoice);
        nt = num2str(length(a));
        if length(a) < mintrialnum %skip if too few trials
            out{i,j,:,:} = [];
        else
            out{i,j,:,:} = plotActivationMap(cPath,animals{i},Rec,{a,b},[animals{i} ' ' Rec ': ' nt ' trials per state'],{'Attentive trials','Bias trials'},clims,fsize,false,false);
            %out{i,j,:,:} = plotVarianceMap(cPath,animals{i},Rec,{a,b},[animals{i} ' ' Rec ': ' nt ' trials per state'],{'Attentive trials','Bias trials'},clims,fsize,false);
            a_rewarded = a_rewarded + sum(SessionData.Rewarded(a));
            b_rewarded = b_rewarded + sum(SessionData.Rewarded(b));
            a_unrewarded = a_unrewarded + sum(~SessionData.Rewarded(a));
            b_unrewarded = b_unrewarded + sum(~SessionData.Rewarded(b));
            a_right = a_right + sum(SessionData.ResponseSide(a) == 1);
            b_right = b_right + sum(SessionData.ResponseSide(b) == 1);
            a_left = a_left + sum(SessionData.ResponseSide(a) == 2);
            b_left = b_left + sum(SessionData.ResponseSide(b) == 2);
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
            loc = loc + 1;
        end

    end
end

attendmean = squeeze(mean(attend,1,'omitnan')); %average over animals/sessions
biasmean = squeeze(mean(bias,1,'omitnan')); %average over animals/sessions
combo = cat(4,attendmean,biasmean);


%% plotting the average
%clims = {[-.01 .01],[-.01 .01]}; %for df/f
clims = {[-.001 .001],[-.001 .001]}; %for variance
pltlegend = {'Engaged trials','Bias trials'};
fsize = 29;
set(gca,'FontSize',fsize)
%plttitle = 'Activity map averaged over sessions';
plttitle = '';
figure('units','normalized','outerposition',[0 0 1 1],'PaperSize',[40 40])
%figure
for i = 1:2
    subplot(3,5,1+(i-1)*5);
    mapImg = imshow(squeeze(combo(1,:,:,i)), clims{1});
    colormap(mapImg.Parent,'colormap_blueblackred'); axis image;
    if i == 1
        title('Baseline','FontSize',fsize);
    end
    set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
    %hcb = colorbar;
    %hcb.Title.String = 'dF/F';
    ylabel(pltlegend{i},'FontSize',fsize);

    subplot(3,5,2+(i-1)*5);
    mapImg = imshow(squeeze(combo(2,:,:,i)), clims{1});
    colormap(mapImg.Parent,'colormap_blueblackred'); axis image;
    if i == 1
        title('Trial Initiation','FontSize',fsize);
    end
    set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
    %hcb = colorbar;
    %hcb.Title.String = 'dF/F';

    subplot(3,5,3+(i-1)*5);
    mapImg = imshow(squeeze(combo(3,:,:,i)), clims{1});
    colormap(mapImg.Parent,'colormap_blueblackred'); axis image;
    if i == 1
        title('Stimulus','FontSize',fsize);
    end
    set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
    %hcb = colorbar;
    %hcb.Title.String = 'dF/F';

    subplot(3,5,4+(i-1)*5);
    mapImg = imshow(squeeze(combo(4,:,:,i)), clims{1});
    colormap(mapImg.Parent,'colormap_blueblackred'); axis image;
    if i == 1
        title('Delay','FontSize',fsize);
    end
    set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
    %hcb = colorbar;
    %hcb.Title.String = 'dF/F';

    subplot(3,5,5+(i-1)*5);
    mapImg = imshow(squeeze(combo(5,:,:,i)), clims{1});
    colormap(mapImg.Parent,'colormap_blueblackred'); axis image;
    if i == 1
        title('Response','FontSize',fsize);
    end
    set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
    hcb = colorbar;
    hcb.Title.String = 'dF/F';
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
hcb.Title.String = 'dF/F Difference';
hcb.Position = hcb.Position + [0.01 0 0 0];
hcb.FontSize = fsize;
sgtitle(plttitle);

%% load the allen map, select regions to plot
clear z zname
t = load('..\lib\allenDorsalMapSM.mat');
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

%% plot PSTHs
for i = 1:length(z)
    A = []; B = [];
    for j = 1:length(animals)
        for k = 1:length(sessiondates{j})
            fprintf('\n%s, %s',animals{j},sessiondates{j}{k})
            [~,inds{1},inds{2}] = getStateInds(cPath,animals{j},sessiondates{j}{k},method,glmFile,dualcase);

            if length(inds{1}) < mintrialnum
                fprintf('\nSkipping!\n');
                continue
            end

            [temp,eventframes] = plotRegionPSTH(cPath,animals{j},sessiondates{j}{k},inds,z(i),zname(i),{'','','','','','','Engaged trials','','','','','','','Biased trials'},false,false);
            exportgraphics(gcf,strjoin(['C:\Data\churchland\PowerpointsPostersPresentations\SFN2022\EMX_individual_psth\' animals{j} sessiondates{j}{k}  zname(i) '.pdf']));
            close gcf;
            A = [A,temp{1}]; B = [B,temp{2}]; %[nframes, ntrials]
        end
    end
    figure;hold on; title(zname(i))
    time = (0:1:size(A,1)-1) ./ 30; %IMPORTANT, FS IS 15 FOR CSTR MICE
    stdshade(A',.2,'red',time,6,eventframes,[]);
    stdshade(B',.2,'blue',time,6,eventframes,[]);
    legend('','','','','','Engaged trials','','','','','','Biased trials');
    xlabel('Time (s)')
    ylabel('dF/F')
    set(gca,'TickDir','out')
    
end



