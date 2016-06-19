% modified BFGS update, refer to 'P.K.H. Phua and R. Setiono. Combined quasi-newton updates for unconstrained optimization'
% vk : f_{k+1} - f_k
% s_k: x_{k+1} - x_k
% y_k: \nabla f_{k+1} - \nabla f_{k}
% H_k: Hessian approximation
% gradplus: \nalba f_{k+1}
%%
function [H ok] = BIGGS(vk, sk, yk, gradplus, Hk)
    ok = true;
    snorm = norm(sk);
    ynorm = norm(yk);
    temp = sk' * yk;
    if temp > sqrt (OptimizeConstant.DBL_EPSILON) * snorm * ynorm
        tk = 2*(sk'*gradplus - vk)/temp;
        H = Hk + 1/temp*((1/tk+(yk'*Hk*yk)/temp) * (sk * sk') - sk * yk' * Hk - Hk * yk *sk');
    else 
        ok = false;
        H = Hk;
    end
end