% Get e^V(k)
%%
function [expV, boolV0] = getExpV(M)
    boolV0 = 1;
    [n m]= size(M);
    b = zeros(n,1);
    b(n)=1;
    b_sp = sparse(b);    
    I = speye(size(M));  
    A = I - M;   
    Z = A\b_sp;
    % Z = linsolve(A,b_sp)
    minele = min(Z);
%     if minele == 0 || minele < OptimizeConstant.NUM_ERROR
%         boolV0 = 0;
%     end    
    if minele < -1e-10 || minele < OptimizeConstant.NUM_ERROR
        boolV0 = 0;
    end    

    Zabs = abs(Z); % MAYBE SET TO VERY SMALL POSITIVE VALUE?  
    resNorm = norm(A * Zabs - b_sp);
    if resNorm > OptimizeConstant.RESIDUAL
        boolV0 = 0;
    end   
    expV = full(Zabs);
end
