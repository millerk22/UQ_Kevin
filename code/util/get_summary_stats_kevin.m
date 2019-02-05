function summary_stats = get_summary_stats_kevin(fid, class_result, ground_truth, ratio)
    fid_list = get_fid_list(fid);
    to_eval = setdiff(ratio+1:numel(ground_truth), fid_list); % get the "test" indices that were not accidentally included in the fid_list?
    assert(numel(to_eval) == numel(ratio+1: numel(ground_truth))); % if to_eval isn't the same number of elements as the testing raise an error
    summary_stats.C = Nind2vec(ground_truth(to_eval))' * Nind2vec(class_result(to_eval));
    summary_stats.recall_confmat = confusion_matrix(ground_truth(to_eval), Nind2vec(class_result(to_eval)));
    summary_stats.precision_confmat = confusion_matrix(class_result(to_eval),Nind2vec(ground_truth(to_eval))); 
    summary_stats.recall = mean(diag(summary_stats.recall_confmat));
    summary_stats.acc = sum(ground_truth(to_eval) == class_result(to_eval))/numel(ground_truth(to_eval:end));
    summary_stats.precision = mean(diag(summary_stats.precision_confmat));
end
