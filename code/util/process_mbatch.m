% process mbatch results
function [C_recall, C_prec] = process_mbatch(filename,batch_size)
    load([filename,'.mat']);
    %%
    N = batch_size * numel(temp);
    result = [];
    new_gt = new_gt(1:N);
    for i = 1:numel(temp)
        result = [result;temp{i}];
    end
    C_recall = confusion_matrix(new_gt, Nind2vec(result));
    acc = sum(result == new_gt)/numel(new_gt);
    nystrom_opt
    disp(['Mean recall: ', num2str(mean(diag(C_recall)))]);
    disp(['Accuracy: ', num2str(acc)]);
    C_prec = confusion_matrix(result, Nind2vec(new_gt));
    disp(['Mean precision: ', num2str(mean(diag(C_prec)))]);
end
