%   BFGS approximation
%%
function [H ok] = BFGS(sk, yk, Hk)
    snorm = norm(sk);
    ynorm = norm(yk);
    temp = sk' * yk;
    if temp > sqrt (OptimizeConstant.DBL_EPSILON) * snorm * ynorm        
        H = yk*yk'/(temp) - ((Hk * sk) * (sk' *Hk))/(sk' * Hk * sk) + Hk;
        ok = true;
    else 
        H = Hk;
        ok = false;
    end
    
end