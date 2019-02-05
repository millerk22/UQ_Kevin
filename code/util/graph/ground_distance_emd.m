
%% Compute ground distance matrix
matrix = zeros(8,8);
interval = 2.0 * pi / 4.0;
opposite = [5,6,7,8,1,2,3,4];
sequence = 1 : 8;

for i = 1:8
    opp = opposite(i);
    matrix(i,opp) = 4.0 * interval;
    
    for j = 1 : 3
        if opp - j <= 0
            left = 8 + (opp-j);
        else
            left = opp-j;
        end
        
        if opp + j > 8
            right = opp + j - 8;
        else
            right = opp + j;
        end

        matrix(i,left) = (4-j) * interval;
        matrix(i,right) = (4-j) * interval;
    end
end

% Convert the distance matrix into a vector
coeff_vec = matrix(:);

save('ground_distance.mat','coeff_vec')


