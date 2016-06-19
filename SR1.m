% Rank 1 approximation
%%
function [H ok] = SR1(sk, yk, Hk)
    pk = yk - Hk * sk;
    uk = sk' * pk;
    if abs(uk) > OptimizeConstant.R_SR1 * norm(sk) * norm (pk)
        H = Hk + (pk * pk')/uk ; 
        ok = true;
    else
        H = [];
        ok = false;        
    end
end