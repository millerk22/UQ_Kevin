function ratio =  huji_sampling_60seg(rng_seed)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Roughly 50-50 train-test split
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('./data/HUJI/ground_truth_huji_disney_60seg.mat')
load('./data/HUJI/index_huji_disney_seg.mat')
load('./data/HUJI/H_huji_disney_60seg.mat')
rng(rng_seed);
min_label_per_class = [900,900,1300,1300,1300,1300,1300] * 1;
label = {4,7,2,5,3,6,1};
MIN_PER_CLASS = containers.Map(label,min_label_per_class);

used_list = [];

class_count = zeros(7,1);

status = true;

train_index = [];
train_id_list = [];
failed_list = [];
nVid = numel(index);
u = 0;

for k = 1 : 7
    class = label{k};
    
    failed_list = [];
    
    while class_count(class) < MIN_PER_CLASS(class)
        not_train_yet = setdiff(1:nVid,train_id_list);
        
        if numel(failed_list) == numel(not_train_yet)
            fprintf('Not enough training dataset for class %d\n',class)
            break
        end
        not_train_yet = setdiff(not_train_yet,failed_list);
        choosen_seq_ind = not_train_yet(randi(numel(not_train_yet),1,1));
        start = index(choosen_seq_ind)+1;
        
        if choosen_seq_ind == nVid
            finish = numel(ground_truth);
        else
            finish = index(choosen_seq_ind+1);
        end
        
        gt = ground_truth(start:finish);
        
        
        count = numel(find(gt==class));
        
        if  count == 0
            failed_list = [failed_list,choosen_seq_ind];
        else
            train_id_list = [train_id_list,choosen_seq_ind];
            train_index = [train_index,start:finish];
            class_count(class) = class_count(class) + count;
        end
         
    end
    
end

    

    

%disp(class_count)
%disp(sum(class_count))



test_index = setdiff(1:size(ground_truth,1),train_index);

train_gt = ground_truth(train_index);
train_H = H(:,train_index);

test_gt = ground_truth(test_index);
test_H = H(:,test_index);

H = [train_H,test_H];
ground_truth = [train_gt;test_gt];
ratio = numel(train_gt);%/numel(ground_truth);

%save('ground_truth_huji_sampling_60seg.mat','ground_truth','-mat')
%save('H_huji_sampling_60seg.mat','H','-mat')

save('./data/HUJI//temp/ground_truth.mat','ground_truth','ratio','-mat')
save('./data/HUJI/temp/H.mat','H','-mat')
save('./data/HUJI/temp/train_test_index.mat', 'train_index','test_index')