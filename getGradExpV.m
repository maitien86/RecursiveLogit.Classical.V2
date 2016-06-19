function [gradExpV] = getGradExpV(M, Att, op, expV)
    I = speye(size(M));  
    A = I - M; 
    gradExpV = zeros(size(M,1),op.n);
    for i = 1:op.n
        u = M .* (Att(i).Value); 
        v = u * expV; 
        gradExpV(:,i) = A\v; 
    end
end
