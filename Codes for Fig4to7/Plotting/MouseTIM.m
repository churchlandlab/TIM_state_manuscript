mouse = 'JC047';


D = load(['X:\StateProjectCentralRepo\JoaoDataset_GRB\', mouse, '\saved_data\', mouse, '_analyzed_data.mat']);
rawDLC = load(['X:\StateProjectCentralRepo\JoaoDataset_GRB\', mouse, '_raw_DLC_data.mat']);
rawPupil = load(['X:\StateProjectCentralRepo\JoaoDataset_GRB\', mouse, '_raw_pupil.mat']);
rawPupil = rawPupil.raw_pupil;

num_Session = size(rawPupil, 1);





%% Linear regression for TIM calculatiuon 
%% We only use the period between stimulus on and spouts coming in
%% The camera frame rates can be 30Hz or 60Hz

all_Sessions = {};
all_Pupil = {}; 

for i = 1 : num_Session
    
    frameRate = 3000 / (D.aligned_frametimes{i, 1}.frame_times(3000) - D.aligned_frametimes{i, 1}.frame_times(1));
    sessionLength = size(D.eventFrames{i,1}.trialstart_frames, 1);
    
    stimOn_Frame = arrayfun(@(x)(find(D.aligned_frametimes{i,1}.frame_times >= x, 1)), ...
        D.aligned_frametimes{i,1}.stim_onsets);
    
    
    % "response_times" provides the time difference between "trialstart_times" and animal's 2nd lick.
    % 1) there is a pre-stimulus delay between "trialstart_times" and "stim_onsets"
    % 2) in some trials, licking may start slightly before the end of stimulus, we here assume this gap is 200ms long at maximum
    stimulus_time = D.aligned_frametimes{i,1}.response_times - (D.aligned_frametimes{i,1}.stim_onsets - D.aligned_frametimes{i,1}.trialstart_times);
    stimulus_time = stimulus_time - 0.2;
    
    stimulus_frame = floor(stimulus_time .* frameRate);
    
    
    % the time length of stimulus is variable (1.0~1.5s). We only use the 70% length of the longest trial
    a = sort(stimulus_frame);
    a = a(~isnan(a));
    a = a(round(length(a) * 0.7));
    
    DLC_L = rawDLC.DLC_Lateral_raw{i, 1};
    DLC_L(:, 3:3:end) = [];
    DLC_B = rawDLC.DLC_Bottom_raw{i, 1};
    DLC_B(:, 1:6) = [];       % The first two labels in the bottom camera are water spouts.
    DLC_B(:, end-2:end) = []; % The last label is tongue, which has very low confidence because it isn't shown in most of the frames.
    DLC_B(:, 3:3:end) = [];
    
    
    thisSession_L = nan(sessionLength, a, size(DLC_L, 2));
    thisSession_B = nan(sessionLength, a, size(DLC_B, 2));
    
    for ii = 1 : sessionLength
        stimulusEnd = stimulus_frame(ii);
        if stimulusEnd >= a 
            thisSession_L(ii, :, :) = DLC_L(stimOn_Frame(ii) : stimOn_Frame(ii)+a-1, :);
            thisSession_B(ii, :, :) = DLC_B(stimOn_Frame(ii) : stimOn_Frame(ii)+a-1, :);
        elseif stimulusEnd < a 
            thisSession_L(ii, 1:stimulusEnd, :) = DLC_L(stimOn_Frame(ii) : stimOn_Frame(ii)+stimulusEnd-1, :);
            thisSession_B(ii, 1:stimulusEnd, :) = DLC_B(stimOn_Frame(ii) : stimOn_Frame(ii)+stimulusEnd-1, :);
        end   
    end
    
    
    
    
    
    thisSession = cat(3, thisSession_L, thisSession_B);
    thisSession = thisSession - nanmean(thisSession, 1);
    
    all_Sessions{i} = thisSession;
    
    
    
    ttt = round(frameRate * 0.5);
    thisPupil = rawPupil{i, 1};
    
    all_Pupil{i} = arrayfun(@(x)(nanmean(thisPupil(x - ttt : x - 1))), stimOn_Frame);
   
end




num_bin = cellfun(@(X)(size(X,2)), all_Sessions);
num_bin = min(num_bin);
all_Sessions = cellfun(@(X)(X(:, 1:num_bin, :)), all_Sessions, 'UniformOutput', 0);



num_trials = cellfun(@(X)(size(X,1)), all_Sessions);
DLC_matrix = cat(1, all_Sessions{:});


% Some animals have 60Hz frame rate. We downsample them to 30Hz.
if round(frameRate) == 60
    DLC_matrix = DLC_matrix(:, 1:2:end, :);
    num_bin = size(DLC_matrix, 2);
end


DLC_matrix = areaNormalization(DLC_matrix);
DLC_matrix_size = size(DLC_matrix);

% 5-sigma bar for standardized DLC output
DLC_matrix = reshape(DLC_matrix, [], DLC_matrix_size(3));
for i = 1 : 2 : DLC_matrix_size(3)
    x_std = nanstd(DLC_matrix(:, i));
    y_std = nanstd(DLC_matrix(:, i+1));
    
    x_mean = nanmean(DLC_matrix(:, i));
    y_mean = nanmean(DLC_matrix(:, i+1));
    
    
    idx = [find(abs(DLC_matrix(:, i) - x_mean) > 5*x_std); find(abs(DLC_matrix(:, i+1) - y_mean) > 5*y_std)];
    idx = unique(idx);
    
    DLC_matrix(idx, i) = NaN;
    DLC_matrix(idx, i+1) = NaN;
    
end


DLC_matrix = reshape(DLC_matrix, DLC_matrix_size);





c = cumsum(num_trials);
c = [0, c];
for i = 1 : num_Session 
    all_Sessions{1, i} = DLC_matrix(c(i)+1:c(i+1), :, :);   
end




X = D.group_info;
a = ((X(:, 3) == X(:, 4)) - 0.5) * 2;
a(X(:, 4) == 0) = 0;

b = ((X(:, 5) == X(:, 6)) - 0.5) * 2;
a(X(:, 6) == 0) = 0;

X(X == 1) = -1;
X(X == 2) = 1;

X = [X, a, b];
X(:, 1) = [];   % No opto trials in this dataset

XX = {};
for i = 1 : num_Session 
    XX{1, i} = X(c(i)+1:c(i+1), :);   
end
X = XX;
clear XX



Predict = {};
for i = 1 : num_Session
    t = nan(size(all_Sessions{1, i}));
    
    for ii = 1 : num_bin
        
        for iii = 1 : size(DLC_matrix, 3)
            y = all_Sessions{1, i}(:, ii, iii);
            y_idx = find(~isnan(y));
            
            mdl = fitrlinear(X{1, i}, y, 'Kfold', 5);
            t(y_idx, ii, iii) = kfoldPredict(mdl);
                
        end
        
    end
    Predict = [Predict, t];
end


Predict_backup = Predict;
Predict = cat(1, Predict{:});






RSS = (Predict - DLC_matrix) .^ 2;
RSS = squeeze(nanmean(RSS, 1));

TSS = (DLC_matrix - nanmean(DLC_matrix, 1)) .^ 2;
TSS = squeeze(nanmean(TSS, 1));
R_squared = ones(size(RSS)) - RSS ./ TSS;
R_squared = R_squared';



figure;
hold on
for i = 1 : size(DLC_matrix, 3)
    plot(R_squared(i, :));
end
ylabel('R_squared');
ylim([0 1]);
title(['The Averaged R squared of all labels, ', mouse]);
xlabel('Frame Number');
set(gca,'box','off');
set(gca,'tickdir','out');
hold off






%% TIM & pupil analyses
smooth_window = 50;


CorrectRate = D.trialData.rewarded;
CorrectRate = smoothdata(CorrectRate, 'gaussian', smooth_window);





DLCEnergy = diff(DLC_matrix, 1, 2);
DLCEnergy = sqrt(DLCEnergy(:, :, 1:2:end) .^ 2 + DLCEnergy(:, :, 2:2:end) .^ 2);

B = permute(DLCEnergy,[3 2 1]); 
B_size = size(B);
B = reshape(B, B_size(1), []);

B = normalize(B, 2);
B = reshape(B, B_size);

DLCEnergy = permute(B,[3 2 1]);
clear B B_size

DLCEnergy_backup = DLCEnergy;
% DLCEnergy = DLCEnergy(:, :, [1:5, 10, 15:21, 23:24]);

DLCEnergy = nanmean(DLCEnergy, 3);
DLCEnergy = nanmean(DLCEnergy, 2);
DLCEnergy = smoothdata(DLCEnergy, 'gaussian', smooth_window);







Diff_toRegression = DLC_matrix - Predict;
Diff_toRegression = sqrt(Diff_toRegression(:, :, 1:2:end) .^ 2 + Diff_toRegression(:, :, 2:2:end) .^ 2);

B = permute(Diff_toRegression,[3 2 1]); 
B_size = size(B);
B = reshape(B, B_size(1), []);

B = normalize(B, 2);
B = reshape(B, B_size);

Diff_toRegression = permute(B,[3 2 1]);
clear B B_size



TIM = nanmean(Diff_toRegression,3);
TIM = nanmean(TIM,2);

TIM = smoothdata(TIM, 'gaussian', smooth_window);





figure('name', 'Performance, TIM, and Motion Energy');
hold on
plot(TIM);
plot(DLCEnergy);
ylabel('Variance/Motion Energy');

yyaxis right
plot(CorrectRate');
ylabel('Correct Rate');
ylim([0.5 1.0]);
xlabel('Trial Number');



if length(c) > 1
    for i = 1 : (length(c)-1) % plotting session borders
        line([c(i) c(i)], [0.5 1], 'Color','black','LineStyle','--');
    end
end
legend({'TIM','DLC Energy','Performance'});
set(gca,'box','off');
set(gca,'tickdir','out');

hold off






Pupil = cat(1, all_Pupil{:});
Pupil = smoothdata(Pupil, 'gaussian', smooth_window);


CorreCoeffs = [TIM, DLCEnergy, CorrectRate, Pupil];
CorreCoeffs = corrcoef(CorreCoeffs, 'Rows','complete');
for i = 1 : size(CorreCoeffs, 1)
    CorreCoeffs(i, i) = NaN;
end
figure('Name', 'All the CorreCoeffs');
h = heatmap(CorreCoeffs);
h.XDisplayLabels = {'TIM', 'DLC Energy (normalized)','Performance', 'Pupil'};
h.YDisplayLabels = {'TIM', 'DLC Energy (normalized)','Performance', 'Pupil'};




%%
Results.Pupil = Pupil;
Results.TIM = TIM;
Results.DLCEnergy = DLCEnergy;
Results.Performance = CorrectRate;
Results.rewarded = D.trialData.rewarded;
Results.Pupil_raw = all_Pupil;
Results.DLC = all_Sessions;
Results.Prediction = Predict_backup;






