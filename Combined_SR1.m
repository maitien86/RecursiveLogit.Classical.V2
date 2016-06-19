% Combined approximation for Hessian matrix
% with rank 1 update.
%%
function [H ok] = Combined_SR1(sk, yk, Ak)
    global Op;
    Bh = BHHH();
    uk = yk - Bh*sk;
    [Op.Ak ,ok] = SR1(sk, uk, Ak);
    H = Bh + Op.Ak;
end