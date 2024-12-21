function out = returnVarianceMovie(cpath,mouse,rec,modelfile)
datapath = [cpath filesep mouse filesep 'SpatialDisc' filesep rec filesep];
load('C:\Data\churchland\ridgeModel\widefield\allenDorsalMapSM.mat','dorsalMaps'); %Get allen atlas
try
    load([datapath modelfile],'cMovie');
catch
    fprintf('\nThe encoding model file does not exist!');
    out = NaN;
return
end

out = double(cMovie);
out = mean(out,1);

end
