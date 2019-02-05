function summary_stats = get_summary_stats(fid, class_result, ground_truth, ratio)
    fid_list = get_fid_list(fid);
    to_eval = setdiff(ratio+1:numel(ground_truth), fid_list);
    assert(numel(to_eval) == numel(ratio+1: numel(ground_truth)));
    summary_stats.C = Nind2vec(ground_truth(to_eval))' * Nind2vec(class_result(to_eval));
    summary_stats.recall_confmat = confusion_matrix(ground_truth(to_eval), Nind2vec(class_result(to_eval)));
    summary_stats.precision_confmat = confusion_matrix(class_result(to_eval),Nind2vec(ground_truth(to_eval)));
    summary_stats.recall = mean(diag(summary_stats.recall_confmat));
    summary_stats.acc = sum(ground_truth(to_eval) == class_result(to_eval))/numel(ground_truth(to_eval:end));
    summary_stats.precision = mean(diag(summary_stats.precision_confmat));
end
