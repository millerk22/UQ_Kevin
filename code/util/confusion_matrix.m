function [ C ] = confusion_matrix(ground_truth, yp )
% generate confusion matrix 
NClass = size(yp, 2); 
N = size(yp, 1); 
yg = Nind2vec(ground_truth, NClass);
C = yg'*yp / N; 
for i = 1:NClass
    C(i, :) = C(i, :) / sum(C(i, :)); 
end
end

