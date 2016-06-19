% Combined approximation for Hessian matrix
% with BFGS
function [H ok] = Combined_BFGS(sk, yk, Ak)
    global Op;
    Bh = BHHH();
    uk = yk - Bh*sk;
    [Op.Ak ,ok] = BFGS(sk, uk, Ak);
    H = Bh + Op.Ak;
end