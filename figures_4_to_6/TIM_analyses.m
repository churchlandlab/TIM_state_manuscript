
function [CorrectRate_unsmoothed, DLC_matrix, Fitted, DLCEnergy_raw, Distance_TIM_raw, Results, State_results, borders] = TIM_analyses...
    (Lateral_allFrames, Bottom_allFrames, aligned_FrameTime, raw_data, smooth_window, outlier_filter, HMM_states)

%% 1st, showing the distribution of the timelength of different task epochs. (baseline is always 15 frames.)
global mouse_name Laterl_labels Bottom_labels
mousename = mouse_name; 
DLC_labels = [Laterl_labels, Bottom_labels];

length_stim_L = aligned_FrameTime.stimOff(1,:) - aligned_FrameTime.stimOn(1,:);
length_stim_B = aligned_FrameTime.stimOff(2,:) - aligned_FrameTime.stimOn(2,:);

figure('Name', 'Epoch time length distribution');
subplot(2,2,1);
hold on
histogram(length_stim_L);
histogram(length_stim_B);
title('Stimulus');
xlabel('Frame Number');
ylabel('Trial Number');
legend({'Lateral', 'Bottom'});
disp(['The biggest mis-match between lateral and bottom cameras for the stimulus window is ', num2str(max(abs(length_stim_L - length_stim_B)))]);



length_delay_L = aligned_FrameTime.spoutsIn(1,:) - aligned_FrameTime.stimOff(1,:);
length_delay_B = aligned_FrameTime.spoutsIn(2,:) - aligned_FrameTime.stimOff(2,:);

subplot(2,2,2);
hold on
histogram(length_delay_L);
histogram(length_delay_B);
title('Delay');
xlabel('Frame Number');
ylabel('Trial Number');
disp(['The biggest mis-match between lateral and bottom cameras for the delay window is ', num2str(max(abs(length_delay_L - length_delay_B)))]);



length_response_L = cellfun(@(x)(size(x,1)),Lateral_allFrames.frames,'UniformOutput',true) - aligned_FrameTime.spoutsIn(1,:);
length_response_B = cellfun(@(x)(size(x,1)),Bottom_allFrames.frames,'UniformOutput',true) - aligned_FrameTime.spoutsIn(2,:);

subplot(2,2,3);
hold on
histogram(length_response_L);
histogram(length_response_B);
title('Response');
xlabel('Frame Number');
ylabel('Trial Number');
disp(['The biggest mis-match between lateral and bottom cameras for the response window is ', num2str(max(abs(length_response_L - length_response_B)))]);



length_all_L = 15 + length_stim_L + length_delay_L + length_response_L;
length_all_B = 15 + length_stim_B + length_delay_B + length_response_B;

a = [length_all_L,length_all_B];
subplot(2,2,4);
hold on
histogram(length_all_L, [min(a):2:max(a)]);
histogram(length_all_B, [min(a):2:max(a)]);
title('All Epochs');
xlabel('Frame Number');
ylabel('Trial Number');
disp(['The biggest mis-match between lateral and bottom cameras for the entire trial is ', num2str(max(abs(length_all_L - length_all_B)))]);





%% 2nd, setting a bar for epoch chunking, which is the length including 70% trials
% If the cutoff for the two cameras are different, then use the shorter one.
cutoff_stim_L = sort(length_stim_L);
cutoff_stim_L = cutoff_stim_L(floor(length(cutoff_stim_L)*0.7));

cutoff_stim_B = sort(length_stim_B);
cutoff_stim_B = cutoff_stim_B(floor(length(cutoff_stim_B)*0.7));

cutoff_stim = min([cutoff_stim_L, cutoff_stim_B]);



cutoff_delay_L = sort(length_delay_L);
cutoff_delay_L = cutoff_delay_L(floor(length(cutoff_delay_L)*0.7));

cutoff_delay_B = sort(length_delay_B);
cutoff_delay_B = cutoff_delay_B(floor(length(cutoff_delay_B)*0.7));

cutoff_delay = min([cutoff_delay_L, cutoff_delay_B]);


clear cutoff_stim_L cutoff_stim_B cutoff_delay_L cutoff_delay_B 







%% 3rd, epoch alignment
num_label_L = size(Lateral_allFrames.frames{1,1}, 2)/3;
num_label_B = size(Bottom_allFrames.frames{1,1}, 2)/3;
num_trial = length(length_stim_L);



DLCposition_L_baseline = nan(num_trial, 15, num_label_L*2);
DLCposition_B_baseline = nan(num_trial, 15, num_label_B*2);
for i = 1 : num_trial
    L1 = aligned_FrameTime.stimOn(1,i);
    B1 = aligned_FrameTime.stimOn(2,i);
    
    a = Lateral_allFrames.frames{1, i};
    a(:, 3:3:end) = [];
    
    b = Bottom_allFrames.frames{1, i};
    b(:, 3:3:end) = [];
    
    
    DLCposition_L_baseline(i, :, :) = a(L1-14:L1, :);
    DLCposition_B_baseline(i, :, :) = b(B1-14:B1, :);
end
DLC_matrix_baseline = cat(3, DLCposition_L_baseline, DLCposition_B_baseline);
clear i a b L1 B1





DLCposition_L_stim = nan(num_trial, cutoff_stim, num_label_L*2);
DLCposition_B_stim = nan(num_trial, cutoff_stim, num_label_B*2);
    
    
for i = 1 : num_trial
    
    
    L1 = aligned_FrameTime.stimOn(1,i);
    B1 = aligned_FrameTime.stimOn(2,i);
    
    L2 = aligned_FrameTime.stimOff(1,i);
    B2 = aligned_FrameTime.stimOff(2,i);
    

    
    if (L2 - L1) >= cutoff_stim
        
        DLC_thisTrial_L = Lateral_allFrames.frames{1, i};
        DLC_thisTrial_L(:, 3:3:end) = [];
        DLC_thisTrial_L = DLC_thisTrial_L(L1:(L1+cutoff_stim-1), :);
        DLCposition_L_stim(i, :, :) = DLC_thisTrial_L;
    else
        a = Lateral_allFrames.frames{1, i};
        a(:, 3:3:end) = [];
        a = a(L1:(L2-1), :);
        
        DLCposition_L_stim(i, 1:size(a,1), :) = a;
    end
    
    
    if (B2 - B1) >= cutoff_stim
        
        DLC_thisTrial_B = Bottom_allFrames.frames{1, i};
        DLC_thisTrial_B(:, 3:3:end) = [];
        DLC_thisTrial_B = DLC_thisTrial_B(B1:(B1+cutoff_stim-1), :);
        DLCposition_B_stim(i, :, :) = DLC_thisTrial_B;
    else
        a = Bottom_allFrames.frames{1, i};
        a(:, 3:3:end) = [];
        a = a(B1:(B2-1), :);
        
        DLCposition_B_stim(i, 1:size(a,1), :) = a;
    end
end

DLC_matrix_stim = cat(3, DLCposition_L_stim, DLCposition_B_stim);
clear L1 L2 B1 B2 a i 



DLCposition_L_delay = nan(num_trial, cutoff_delay, num_label_L*2);
DLCposition_B_delay = nan(num_trial, cutoff_delay, num_label_B*2);


for i = 1 : num_trial
    
    
    L1 = aligned_FrameTime.stimOff(1,i);
    B1 = aligned_FrameTime.stimOff(2,i);
    
    L2 = aligned_FrameTime.spoutsIn(1,i);
    B2 = aligned_FrameTime.spoutsIn(2,i);
    

    if (L2 - L1) >= cutoff_delay
        
        DLC_thisTrial_L = Lateral_allFrames.frames{1, i};
        DLC_thisTrial_L(:, 3:3:end) = [];
        DLC_thisTrial_L = DLC_thisTrial_L(L1:(L1+cutoff_delay-1), :);
        DLCposition_L_delay(i, :, :) = DLC_thisTrial_L;
    else
        a = Lateral_allFrames.frames{1, i};
        a(:, 3:3:end) = [];
        a = a(L1:(L2-1), :);
        
        DLCposition_L_delay(i, 1:size(a,1), :) = a;
    end
    
    
    if (B2 - B1) >= cutoff_delay
        
        DLC_thisTrial_B = Bottom_allFrames.frames{1, i};
        DLC_thisTrial_B(:, 3:3:end) = [];
        DLC_thisTrial_B = DLC_thisTrial_B(B1:(B1+cutoff_delay-1), :);
        DLCposition_B_delay(i, :, :) = DLC_thisTrial_B;
    else
        a = Bottom_allFrames.frames{1, i};
        a(:, 3:3:end) = [];
        a = a(B1:(B2-1), :);
        
        DLCposition_B_delay(i, 1:size(a,1), :) = a;
    end
end

DLC_matrix_delay = cat(3, DLCposition_L_delay, DLCposition_B_delay);
clear L1 L2 B1 B2 a i




t = min([length_response_L, length_response_B]);

DLCposition_L_response = nan(num_trial, t, num_label_L*2);
DLCposition_B_response = nan(num_trial, t, num_label_B*2);
for i = 1 : num_trial
    L1 = aligned_FrameTime.spoutsIn(1,i);
    B1 = aligned_FrameTime.spoutsIn(2,i);
    
    a = Lateral_allFrames.frames{1, i};
    a(:, 3:3:end) = [];
    
    b = Bottom_allFrames.frames{1, i};
    b(:, 3:3:end) = [];
    
    
    DLCposition_L_response(i, :, :) = a(L1 : L1+t-1, :);
    DLCposition_B_response(i, :, :) = b(B1 : B1+t-1, :);
end
DLC_matrix_response = cat(3, DLCposition_L_response, DLCposition_B_response);
clear i a b L1 B1 t






%% 4th, preparing the loading matrix
DLC_matrix = cat(2, DLC_matrix_stim, DLC_matrix_delay);
% DLC_matrix = DLC_matrix_delay;
% DLC_matrix = areaNormalization(DLC_matrix);



temp = reshape(DLC_matrix, [num_trial, size(DLC_matrix,2)*(num_label_B+num_label_L)*2]);
variableList = {'optotype', 'stimulus', 'outcome', 'responseside', 'formeroutcome', 'formerresponseside'};
for a = 1 : length(variableList)
    temp = GroupData_Mouse(temp, variableList{a}, raw_data);
end

group_info = temp(:, end-length(variableList)+1:end);
clear a temp



factorTime = zeros(num_trial, size(DLC_matrix,2), 7);    %7 task variables


for a = 1 : num_trial
    
    % stimulus
    factorTime(a, :, 1) = (group_info(a, 2) - 1.5) * 2;
    
    
    % outcome 
    factorTime(a, :, 2) = (group_info(a, 3) - 1.5) * 2;

    
    % responseside (taking the choice bias into account, this variable starts at the very begining of the trial)
    factorTime(a, :, 3) = (group_info(a, 4) - 1.5) * 2;
    
    % interaction between responside and outcome (trial n)
    if group_info(a,4) == 1 && group_info(a, 3) == 1
        factorTime(a, :, 4) = 1;    % left choice & correct
    elseif group_info(a,4) == 2 && group_info(a, 3) == 1
        factorTime(a, :, 4) = -1;    % right choice & correct
    elseif group_info(a,4) == 1 && group_info(a, 3) == 2
        factorTime(a, :, 4) = -1;    % left choice & wrong
    elseif group_info(a,4) == 2 && group_info(a, 3) == 2
        factorTime(a, :, 4) = 1;    % right choice & wrong
    end
    
    
    % formeroutcome
    factorTime(a, :, 5) = (group_info(a,5) - 1.5) * 2;

    
    %formerresponseside
    factorTime(a, :, 6) = (group_info(a,6) - 1.5) * 2;
    
    % interaction between responside and outcome (trial n-1)
    if group_info(a,6) == 1 && group_info(a, 5) == 1
        factorTime(a, :, 7) = 1;
    elseif group_info(a,6) == 2 && group_info(a, 5) == 1
        factorTime(a, :, 7) = -1;
    elseif group_info(a,6) == 1 && group_info(a, 5) == 2
        factorTime(a, :, 7) = -1;
    elseif group_info(a,6) == 2 && group_info(a, 5) == 2
        factorTime(a, :, 7) = 1;
    end
    
end

factorTime(factorTime == -3) = 0;   % -3 are non-previous-choice trials
factorTime = squeeze(factorTime(:, 1, :));



% Applying the filters and remove the too short trials (some delay windows are shorter than 4 frames)
% removing no-reponses trials and opto trials in case they are outliers
deleted = logical((group_info(:,4) == 0) + (group_info(:,6) == 0) + (group_info(:,1) ~= 0));

% Remove the outliers if applying the filter
DLC_matrix_size = size(DLC_matrix);
DLC_matrix = reshape(DLC_matrix, [], DLC_matrix_size(3));


if outlier_filter > 0  
    for i = 1 : 2 : DLC_matrix_size(3)
        x_std = nanstd(DLC_matrix(:, i));
        y_std = nanstd(DLC_matrix(:, i+1));
        
        x_mean = nanmean(DLC_matrix(:, i));
        y_mean = nanmean(DLC_matrix(:, i+1));
        
        
        idx = [find(abs(DLC_matrix(:, i) - x_mean) > outlier_filter*x_std); find(abs(DLC_matrix(:, i+1) - y_mean) > outlier_filter*y_std)];
        idx = unique(idx);
        
        DLC_matrix(idx, i) = NaN;
        DLC_matrix(idx, i+1) = NaN;        
    end   
end
DLC_matrix = reshape(DLC_matrix, DLC_matrix_size);
clear i x_std x_mean y_std y_mean


window_length = DLC_matrix(:,:,1);
window_length = isnan(window_length);
window_length = size(window_length, 2) - sum(window_length,2);
% Delete the trials shorter than 30 frames.
deleted(window_length<30) = 1;




factorTime(deleted, :) = [];
DLC_matrix(deleted, :, :) = [];

% Balance the number of engaged/disengaged trials before training the model
if ~isempty(HMM_states)
    HMM_states(deleted) = [];

    if sum(isnan(HMM_states)) > 0
        HMM_states = fillmissing(HMM_states, 'linear');
    end
end





%% 5th, regression
borders = raw_data.nTrials; % for plotting session borders
ttt = [0, cumsum(raw_data.nTrials)];
for i = 1 : length(raw_data.nTrials)
    borders(i) = borders(i) - sum(deleted((ttt(i)+1) : ttt(i+1)));
end
clear ttt i
borders = [1, cumsum(borders)];



Fitted = nan(num_trial - sum(deleted), size(DLC_matrix,2), (num_label_B+num_label_L)*2);

for e = 1 : (length(borders) - 1)
    
    idx = [borders(e) : borders(e+1)];
    idx_correct = idx(factorTime([borders(e) : borders(e+1)], 2) == 1);
    idx_wrong = idx(factorTime([borders(e) : borders(e+1)], 2) == -1);
    
    if length(idx_correct) > length(idx_wrong)
        idx_correct = idx_correct(randsample(length(idx_correct), length(idx_wrong)));
    else
        idx_wrong = idx_wrong(randsample(length(idx_wrong), length(idx_correct)));
    end
    
    
    
    idx_balanced = [idx_correct, idx_wrong];
    
    idx_extra = setdiff([borders(e) : borders(e+1)], idx_balanced);
    

    DLC_training = DLC_matrix(idx_balanced, :, :);
    X = factorTime(idx_balanced, :);
    
    idx_all = [idx_balanced, idx_extra];
    X_all = factorTime(idx_all, :);
    
    for a = 1 : size(DLC_training,2)
        
        thisFrame = squeeze(DLC_training(:, a, :));
        
        
        for b = 1 : (num_label_B+num_label_L)*2
            
            y = thisFrame(:, b);
            if sum(isnan(y)) > 0
                y(isnan(y)) = nanmean(y);
            end
            
            
            mdl = fitlm(X, y);
            Fitted(idx_all, a, b) = predict(mdl, X_all);

            
            clear mdl
            
        end      
    end
end

clear a b c e thisFrame PCA_explained X







%% 6th, TIM analyses
% correct rate
CorrectRate = raw_data.Rewarded;
CorrectRate(deleted) = [];
CorrectRate_unsmoothed = CorrectRate;
CorrectRate = double(CorrectRate);


% motion energy
DLCEnergy = diff(DLC_matrix, 1, 2);
DLCEnergy = sqrt(DLCEnergy(:, :, 1:2:end) .^ 2 + DLCEnergy(:, :, 2:2:end) .^ 2);

B = permute(DLCEnergy,[3 2 1]); 
B_size = size(B);
B = reshape(B, B_size(1), []);

B = normalize(B, 2, 'zscore');
B = reshape(B, B_size);

DLCEnergy = permute(B,[3 2 1]);
clear B B_size

DLCEnergy_backup = DLCEnergy;

DLCEnergy = nanmean(DLCEnergy, 3);
DLCEnergy = nanmean(DLCEnergy, 2);
DLCEnergy_raw = DLCEnergy;


% TIM
Diff_toRegression = DLC_matrix - Fitted;
Diff_toRegression = sqrt(Diff_toRegression(:, :, 1:2:end) .^ 2 + Diff_toRegression(:, :, 2:2:end) .^ 2);

B = permute(Diff_toRegression,[3 2 1]); 
B_size = size(B);
B = reshape(B, B_size(1), []);

B = normalize(B, 2, 'zscore');
B = reshape(B, B_size);

Diff_toRegression = permute(B,[3 2 1]);
clear B B_size


Distance_TIM_backup = Diff_toRegression;
% Diff_toRegression = Diff_toRegression(:, :, [1:5, 10, 15:21, 23:24]);

Distance_TIM = nanmean(Diff_toRegression,3);
Distance_TIM = nanmean(Distance_TIM,2);
Distance_TIM_raw = Distance_TIM;





for i = 1 : (length(borders) - 1)
    CorrectRate(borders(i):borders(i+1)) = smoothdata(CorrectRate(borders(i):borders(i+1)), 'gaussian', smooth_window);
    
    DLCEnergy(borders(i):borders(i+1)) = smoothdata(DLCEnergy(borders(i):borders(i+1)), 'gaussian', smooth_window);
    
    Distance_TIM(borders(i):borders(i+1)) = smoothdata(Distance_TIM(borders(i):borders(i+1)), 'gaussian', smooth_window);
    
    if ~isempty(HMM_states)
        HMM_states(borders(i):borders(i+1)) = smoothdata(HMM_states(borders(i):borders(i+1)), 'gaussian', smooth_window);
    end
end



State_results = State_TIM_Motionenergy(HMM_states, Distance_TIM, DLCEnergy, borders(2:end), smooth_window, Distance_TIM_raw);
State_results.HMM_state = HMM_states;



end