function [dlc_labels, event_frames] = loadDLCLabels(dpath, animal, rec, camera)
data = load([dpath filesep animal '.mat']);
rec = datetime(rec,'InputFormat','dd-MMM-yyyy');
rec.Format = 'MMMdd_yyyy';

names = fieldnames(data.(animal));
ind = contains(names, string(rec));
name = names{ind};


data = data.(animal).(name);
event_frames = data.aligned_FrameTime;
dlc_labels = data.(camera);

fn = fieldnames(event_frames);
for k=1:numel(fn)
    if( isnumeric(event_frames.(fn{k})) )
        if strcmp(camera,'Lateral')
            temp = event_frames.(fn{k});
            event_frames.(fn{k}) = temp(1,:);
        elseif strcmp(camera,'Bottom')
            temp = event_frames.(fn{k});
            event_frames.(fn{k}) = temp(2,:);
        end
    end
end
end