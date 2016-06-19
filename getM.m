%   Get MUtility
%%
function Mfull = getM(x,AttsLc)   
    global incidenceFull;
    global Op;
    u = AttsLc(1).value * x(1); 
    for i = 2:Op.n
        u = u + AttsLc(i).value * x(i);
    end
    expM = u;
    expM(find(incidenceFull)) = exp(u(find(incidenceFull)));
    Mfull = incidenceFull .* expM;
end