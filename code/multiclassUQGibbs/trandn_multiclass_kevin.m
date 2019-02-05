classdef trandn_multiclass_kevin
    properties
        K
        Qs
        inv_Qs
    end
    methods
        function obj = trandn_multiclass_kevin(K)
            obj.K = K;
            % With knowing the number of classes, we can already define
            % what our Q matrix will be for the LQ factorization of linear
            % constraints that define our problem v_ai(i) >= v_l(i).
            obj.Qs = {};
            obj.inv_Qs = {};
            for ind = 1 : K
                Q = zeros(K);
                Q(:, 1) = ones(K, 1);
                count = 2;
                for j = 1 : K
                    if j ~= ind
                        Q(ind, count) = 1;
                        Q(j, count) = -1;
                        count = count + 1;
                    end
                end
                disp(ind)
                Q
                disp('  ')
                Q = Q';
                for i = 1 : K
                    Q(i, :) = Q(i, :) / norm(Q(i, :));
                end
                obj.Qs{ind} = Q;
                obj.inv_Qs{ind} = inv(Q);
                disp('    ')
            end
        end
        function v = gen_samples(obj, u, gamma, ind)
            u = u';
            N = size(u, 2);
            Q = obj.Qs{ind};
            Q_inv = obj.inv_Qs{ind};
            u_e = Q * u;        
            v_e = u_e;
            v_e(1, :) = u_e(1, :) + gamma * randn(1, N);
            tmp = u_e(2:end, :);
            tmp = tmp(:);
            u_inf = inf * ones(size(tmp));
            v_tmp = trandn(-tmp/gamma, u_inf)*gamma + tmp;
            v_tmp = reshape(v_tmp, obj.K - 1, N);
            v_e(2:end, :) = v_tmp;
            v = Q_inv * v_e;
            v = v';
        end
    end
end
