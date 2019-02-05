function temp2 = recover_batch_label(gt_data,  temp)
    type = unique(gt_data);
    temp2 = temp;
    for i = 1:numel(type)
        temp2(temp == i) = type(i); 
    end
end