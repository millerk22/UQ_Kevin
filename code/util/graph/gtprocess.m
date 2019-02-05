function ground_truth = gtprocess(raw_gt, seg_len)
%segment the raw groundtruth. 
%each segment contains seg_len frames,the length of the output vector is
%len(raw_gt)/seg_len
%the output groundtruth is the most frequent label of the segment
     gt_len = ceil(numel(raw_gt)/seg_len);%
     ground_truth = zeros(gt_len, 1);
     if gt_len >1
         for ind = 1:(gt_len-1)
             temp = raw_gt((ind-1)*seg_len+1:ind*seg_len);% current segment
             elem = unique(temp);% unique label for a segment
             occur_time = 0;
             for i = 1:numel(elem)       
                 if ( sum((temp(:)==elem(i)))>occur_time)% find the most frequently occurred label
                 occur_time = sum(find(temp(:)==elem(i)));
                 ground_truth(ind) = elem(i);
                 end
             end
         end
     end
     % last segment
             temp = raw_gt((gt_len-1)*seg_len+1:end);% current segment
             elem = unique(temp);% unique label for a segment
             occur_time = 0;
             for i = 1:numel(elem)       
                 if ( sum((temp(:)==elem(i)))>occur_time)% find the most frequently occurred label
                 occur_time = sum(find(temp(:)==elem(i)));
                 ground_truth(end) = elem(i);
                 end
             end
end