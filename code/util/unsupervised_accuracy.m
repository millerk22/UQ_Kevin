function acc = unsupervised_accuracy(ind, ground_truth)
lbs = unique(ground_truth);
N = numel(ground_truth);
NCorrect = 0;
for ind1 = 1:numel(lbs)
    target_class = ground_truth(ind == lbs(ind1));
    lbs2 = unique(target_class);
    Count = 0;
    for ind2 = 1:numel(lbs2)
        Count = max(Count, sum(target_class == lbs2(ind2)));
    end
    NCorrect = NCorrect + Count;
end
acc = NCorrect / N;
end