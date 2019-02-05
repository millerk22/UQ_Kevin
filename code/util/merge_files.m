% merge files
% Combine LAPD_10
H_combined = [];
gt_combined = [];
for i = 1:9
    load(['./data/H_matrix/',num2str(i),'.mat']);
    H_combined = [H_combined H];
    load(['./data/gt_folder/',num2str(i),'.mat']);
    gt_combined = [gt_combined; gt];
end
window_size = 60;
H_combined = movingmean(H_combined',window_size,1)';
clear gt H i 
% convert class labels 
unique_classes = unique(gt_combined);
new_gt = gt_combined;
for i = 1:length(unique_classes)
    new_gt (new_gt == unique_classes(i)) = i; % new_gt 1-13 instead of missing 8 etc.
end
clear i 
save('./data/combined_9_videos.mat');
