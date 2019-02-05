function [ distance_matrix ] = histdist( X,Y )

% vectorized matrix to calculate histogram distance without using for loop

[numX,~] = size(X);
[numY,~] = size(Y);

distance_matrix = zeross(numX,numY);

tic
for i = 1 : numX
    for j = i+1: numY
        W1 = X(i,:);
        W2 = Y(j,:);

        W = [W1;W2];

        [Min,~] = min(W);
        [Max,~] = max(W);

        vec = Min ./ Max;
        vec(isnan(vec)) = 0;
        distance_matrix(i,j) = exp(-sum(ones(size(vec))-vec));
        distance_matrix(j,i) = exp(-sum(ones(size(vec))-vec));
    end
end

toc

save('hist_dist.mat','distance_matrix')








end

