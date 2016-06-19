%   Get Utility
%%
function Ufull = getU(x,AttsLc)
    global incidenceFull;
    global Op;
    u = AttsLc(1).value * x(1); 
    for i = 2:Op.n
        u = u + AttsLc(i).value * x(i);
    end
    Ufull = incidenceFull .* u;
end


