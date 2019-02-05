function [distance_matrix] = emddist(X,Y)

% Reference: http://homepages.inf.ed.ac.uk/rbf/CVonline/LOCAL_COPIES/RUBNER/emd.htm
% Input:
% X,Y: Each row is a feature vector

display('Computing EMD')
m = 50;
n= 50;
[numX,~] = size(X);
[numY,~] = size(Y);
distance_matrix = zeros(numX,numY);

% d = load('ground_distance.mat');
% d = d.coeff_vec;

% Each column of W indicates centroid position
W = load('../../../Scripts/nmf/W.mat');
W = W.W;

[nr,~] = size(W);

d = zeros(m,n);
for i = 1 : m
    for j = 1 : n
        d(i,j) = norm((W(:,i)-W(:,j)), 2);
    end
end


tic

for i = 1 : numX
    for j = i+1 : numY 
        
        [e,~] = emd_mex(X(i,:),Y(j,:),d);
       
        distance_matrix(i,j) = e;
        distance_matrix(j,i) = e;
        
    end
end
        
 t = toc;
save('distance_matrix.mat','distance_matrix')
 display(sprintf('Computing EMD takes %f',t))

end





