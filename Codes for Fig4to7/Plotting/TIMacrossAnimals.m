mSM63 = load('X:\Chaoqun\Papers&Talks\TIV-state Manuscript\Major_TIM_figure\mSM63_No2DNormalization.mat');
mSM64 = load('X:\Chaoqun\Papers&Talks\TIV-state Manuscript\Major_TIM_figure\mSM64_No2DNormalization.mat');
mSM65 = load('X:\Chaoqun\Papers&Talks\TIV-state Manuscript\Major_TIM_figure\mSM65_No2DNormalization.mat');
mSM66 = load('X:\Chaoqun\Papers&Talks\TIV-state Manuscript\Major_TIM_figure\mSM66_No2DNormalization.mat');

%%
%
mSM63 = mSM63.DLC_matrix - mSM63.Fitted;
mSM63 = sqrt(mSM63(:, :, 1:2:end) .^ 2 + mSM63(:, :, 2:2:end) .^ 2);

B = permute(mSM63,[3 2 1]); 
B_size = size(B);
B = reshape(B, B_size(1), []);

B = normalize(B, 2);
B = reshape(B, B_size);

mSM63 = permute(B,[3 2 1]);
clear B B_size
mSM63 = nanmean(mSM63,3);

%
mSM64 = mSM64.DLC_matrix - mSM64.Fitted;
mSM64 = sqrt(mSM64(:, :, 1:2:end) .^ 2 + mSM64(:, :, 2:2:end) .^ 2);

B = permute(mSM64,[3 2 1]); 
B_size = size(B);
B = reshape(B, B_size(1), []);

B = normalize(B, 2);
B = reshape(B, B_size);

mSM64 = permute(B,[3 2 1]);
clear B B_size
mSM64 = nanmean(mSM64,3);


%
mSM65 = mSM65.DLC_matrix - mSM65.Fitted;
mSM65 = sqrt(mSM65(:, :, 1:2:end) .^ 2 + mSM65(:, :, 2:2:end) .^ 2);

B = permute(mSM65,[3 2 1]); 
B_size = size(B);
B = reshape(B, B_size(1), []);

B = normalize(B, 2);
B = reshape(B, B_size);

mSM65 = permute(B,[3 2 1]);
clear B B_size
mSM65 = nanmean(mSM65,3);


%
mSM66 = mSM66.DLC_matrix - mSM66.Fitted;
mSM66 = sqrt(mSM66(:, :, 1:2:end) .^ 2 + mSM66(:, :, 2:2:end) .^ 2);

B = permute(mSM66,[3 2 1]); 
B_size = size(B);
B = reshape(B, B_size(1), []);

B = normalize(B, 2);
B = reshape(B, B_size);

mSM66 = permute(B,[3 2 1]);
clear B B_size
mSM66 = nanmean(mSM66,3);




figure;
subplot(1,2,1);
hold on
a = stdshade(mSM63, 0.4, [0 0.4470 0.7410]);
b = stdshade(mSM64, 0.4, [0.8500 0.3250 0.0980]);
c = stdshade(mSM65, 0.4, [0.9290 0.6940 0.1250]);
d = stdshade(mSM66, 0.4, [0.4940 0.1840 0.5560]);

xlabel('Frame from SpoutsIn');
ylabel('TIM');
title('Averaged TIM, with std shade');

set(gca,'box','off');
set(gca,'tickdir','out');

hold off
clear a b c d


subplot(1,2,2);
hold on
a = StandardErrorShade(mSM63, 0.4, [0 0.4470 0.7410]);
b = StandardErrorShade(mSM64, 0.4, [0.8500 0.3250 0.0980]);
c = StandardErrorShade(mSM65, 0.4, [0.9290 0.6940 0.1250]);
d = StandardErrorShade(mSM66, 0.4, [0.4940 0.1840 0.5560]);


legend([a, b, c, d], {'mSM63', 'mSM64', 'mSM65', 'mSM66'});
xlabel('Frame from SpoutsIn');
ylabel('TIM');
title('Averaged TIM, with standard error shade');

set(gca,'box','off');
set(gca,'tickdir','out');

hold off
clear a b c d













%%
for i = 1 : size(mSM63, 1)
    
    if isnan(mSM63(i, 58))
        a = (59 - find(isnan(mSM63(i, :)), 1)) + 1;
        mSM63(i, :) = circshift(mSM63(i, :), a);
    end
    
end


for i = 1 : size(mSM64, 1)
    
    if isnan(mSM64(i, 57))
        a = (58 - find(isnan(mSM64(i, :)), 1)) + 1;
        mSM64(i, :) = circshift(mSM64(i, :), a);
    end
    
end



for i = 1 : size(mSM65, 1)
    
    if isnan(mSM65(i, 58))
        a = (59 - find(isnan(mSM65(i, :)), 1)) + 1;
        mSM65(i, :) = circshift(mSM65(i, :), a);
    end
    
end

for i = 1 : size(mSM66, 1)
    
    if isnan(mSM66(i, 58))
        a = (59 - find(isnan(mSM66(i, :)), 1)) + 1;
        mSM66(i, :) = circshift(mSM66(i, :), a);
    end
    
end





figure;
subplot(1,2,1);
hold on
a = stdshade(mSM63, 0.4, [0 0.4470 0.7410]);
b = stdshade(mSM64, 0.4, [0.8500 0.3250 0.0980]);
c = stdshade(mSM65, 0.4, [0.9290 0.6940 0.1250]);
d = stdshade(mSM66, 0.4, [0.4940 0.1840 0.5560]);

xlabel('Frame from SpoutsIn');
ylabel('TIM');
title('Averaged TIM, with std shade');

set(gca,'box','off');
set(gca,'tickdir','out');

hold off
clear a b c d


subplot(1,2,2);
hold on
a = StandardErrorShade(mSM63, 0.4, [0 0.4470 0.7410]);
b = StandardErrorShade(mSM64, 0.4, [0.8500 0.3250 0.0980]);
c = StandardErrorShade(mSM65, 0.4, [0.9290 0.6940 0.1250]);
d = StandardErrorShade(mSM66, 0.4, [0.4940 0.1840 0.5560]);


legend([a, b, c, d], {'mSM63', 'mSM64', 'mSM65', 'mSM66'});
xlabel('Frame from SpoutsIn');
ylabel('TIM');
title('Averaged TIM, with standard error shade');

set(gca,'box','off');
set(gca,'tickdir','out');

hold off
clear a b c d





