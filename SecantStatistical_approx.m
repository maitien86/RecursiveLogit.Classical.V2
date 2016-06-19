% Secant Statistical approximation for Hessian update.
%% 
function [H ok] = SecantStatistical_approx(sk, yk, secantType)
    Bh = BHHH();
    if secantType == OptimizeConstant.SSA_BFGS
        [H ok] = BFGS(sk,yk,Bh);
    else
        [H ok] = SR1(sk,yk,Bh); 
    end
end