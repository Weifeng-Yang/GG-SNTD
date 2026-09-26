%% 同样，L21范数缩放后，正元素依然为正元素或者0，负元素依然为负元素或者0，
function x = PROXL21(U, lambda,tao)
% 向量化实现版本，更高效
% 求解：min_{Z >= 0} (1/2)||Z - X||_F^2 + λ∑_i ||Z(i,:)||_2
lambda=tao*lambda;
% [m, n] = size(U);

% 应用非负约束
U(U<0)=0;

% 首先对每行应用软阈值（不考虑非负）
row_norms = sqrt(sum(U.^2, 2));
scaling_factors = max(1 - lambda ./ row_norms, 0);
U = scaling_factors .* U;




% 重新计算行范数（仅考虑正分量）
% for i = 1:m
%     pos_idx = Z(i, :) > 0;
%     if any(pos_idx)
%         % 计算正分量的L2范数
%         norm_pos = norm(Z(i, pos_idx), 2);
%         if norm_pos <= lambda
%             Z(i, :) = 0;
%         else
%             % 重新缩放正分量
%             Z(i, pos_idx) = (1 - lambda / norm_pos) * Z(i, pos_idx);
%         end
%     end
% end
x=U;
end

