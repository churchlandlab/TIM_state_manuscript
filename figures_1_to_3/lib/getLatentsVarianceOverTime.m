function [all_latent_vars,eventframes,time] =  getLatentsVarianceOverTime(cPath,Animal,Rec,inds)
for i = 1:length(inds) %will iterate thru each trial set and calculate a PSTH for those
    [alVc,bhv] = align2behavior(cPath,Animal,Rec,inds{i}, true); %pass trialinds based on the sessiondata file, this function will work out the imaging data
    %movie dims are [pixels, frames, trials]

    latents = alVc.all;
    latents_var = var(latents,[],[2,3],'omitnan');
    %latents_var = var(latents,[],2,'omitnan');
    all_latent_vars{i} = squeeze(latents_var);
end
eventframes = alVc.segFrames;
time = (0:1:size(all_latent_vars{1},2)-1) ./ alVc.fs;

end