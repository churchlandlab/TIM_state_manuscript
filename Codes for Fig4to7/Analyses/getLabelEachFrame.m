function [Lateral_allFrames, Bottom_allFrames] = getLabelEachFrame(fPath, mousename)

lateral_path = [fPath 'Lateral\'];
bottom_path = [fPath 'Bottom\'];

lateral_videos = dir([lateral_path mousename '_SpatialDisc_*.csv']);
num_lateral_videos = size(lateral_videos,1);

bottom_videos = dir([bottom_path mousename '_SpatialDisc_*.csv']);
num_bottom_videos = size(bottom_videos,1);

if num_lateral_videos ~= num_bottom_videos
    fprintf(2, 'Lateral and bottom cameras have different video numbers!! \n');
    return;
end


Lateral_allFrames = {};
Bottom_allFrames = {};
for a = 1 : num_lateral_videos
    Lateral_allFrames.frames{a} = ImportDLC_Lateral([lateral_videos(a).folder '\' lateral_videos(a).name], 4, inf);
    Bottom_allFrames.frames{a} = ImportDLC_Bottom([bottom_videos(a).folder '\' bottom_videos(a).name], 4, inf);
    Lateral_allFrames.filename{a} = lateral_videos(a).name;
    Bottom_allFrames.filename{a} = bottom_videos(a).name;
end


end