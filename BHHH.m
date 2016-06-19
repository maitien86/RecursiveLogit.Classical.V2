%   Statistical approximation
%%
function H = BHHH()
    global Gradient;
    global nbobs;
    global Op;
    H = zeros(Op.n); 
    for i = 1:nbobs
        H = H + Gradient(i,:)'*Gradient(i,:);
    end
    H = H/nbobs;
end