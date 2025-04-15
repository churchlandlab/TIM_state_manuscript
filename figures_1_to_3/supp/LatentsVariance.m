clc;clear all;close all;
addpath('..\lib\')
%%
cPath = 'X:\Widefield'; animals = {'mSM63','mSM64','mSM65','mSM66'};
glmFile = 'X:\Widefield\glm_hmm_models\map_all_subjects_targrate.mat';
sessiondates = getGlobalGLMHMMSessions(glmFile); %get sessions with GLM-HMM data
sessiondates{3} = sessiondates{3}(1:11); % drop mSM65 16-Jul-2018 because it could not be aligned to allen atlas

method = 'cutoff';
mintrialnum = 20;
dualcase = 'reward';
onlyLeftChoice = 0;
%dualcase = 'choice_equal';
fsize = 29;
clims = {[-.0003 .0003],[-.0003 .0003]}; %for variance

%% get variance over time of the latents (SVD temporal components)
A_latents = []; B_latents = [];
count=1;
for j = 1:length(animals)
    for k = 1:length(sessiondates{j})
        fprintf('\n%s, %s',animals{j},sessiondates{j}{k})
        [~,inds{1},inds{2}] = getStateInds(cPath,animals{j},sessiondates{j}{k},method,glmFile,dualcase, onlyLeftChoice);

        if length(inds{1}) < mintrialnum
            fprintf('\nSkipping!\n');
            continue
        end

        [all_latent_vars,eventframes,time] = getLatentsVariance(cPath,animals{j},sessiondates{j}{k},inds);
        A_latents(count,:,:) = all_latent_vars{1};
        B_latents(count,:,:) = all_latent_vars{2};
        count = count+1;
    end
end
%% plot
latent_inds = [1,2,3,4,5,6,7,8];

for i=1:length(latent_inds)
    latent_ind = latent_inds(i);
    figure;hold on;
    stdshade(squeeze(A_latents(:,latent_ind,:)),.2,'red',time,6,eventframes,[]);
    stdshade(squeeze(B_latents(:,latent_ind,:)),.2,'blue',time,6,eventframes,[]);
    legend('','','','','','Engaged trials','','','','','','Biased trials');
    xlabel('Time (s)')
    ylabel(['Variance of SVD dimension ',num2str(latent_ind)]);
    set(gca,'TickDir','out');
end

%% plot a histogram for one mouse
[~,inds{1},inds{2}] = getStateInds(cPath,animals{1},sessiondates{1}{1},method,glmFile,dualcase, onlyLeftChoice);
[all_latent_vars,eventframes,time] = getLatentsVariance(cPath,animals{1},sessiondates{1}{1},inds);

