function alignedMap = returnVarianceMap(cpath,mouse,rec,modelfile)
datapath = [cpath filesep mouse filesep 'SpatialDisc' filesep rec filesep];
try
    load([datapath modelfile],'cMap','U');
catch
    fprintf('\nThe encoding model file does not exist!');
    alignedMap = NaN;
    return
end

%out = double(cMap);
try
    load([datapath 'opts3.mat']);
    opts = opts3;clear opts3;
catch
    try
        load([datapath 'opts2.mat']);
    catch
        load([datapath 'opts.mat']);
    end
end

mask = squeeze(isnan(U(:,:,1)));
map = arrayShrink(cMap,mask,'split');

load('C:\Data\churchland\ridgeModel\allenDorsalMapMM.mat');
allenMask = dorsalMaps.allenMask;

alignedMap = alignAllenTransIm(double(map),opts.transParams); %align to allen atlas
alignedMap = alignedMap(:, 1:size(allenMask,2),:);
edgemap = dorsalMaps.edgeMapScaledMax(1:size(alignedMap,1),:); %allen edge map
alignedMap(allenMask == 1) = NaN;
alignedMap(edgemap == 1) = NaN; %apply allen edge map for visualization
end
