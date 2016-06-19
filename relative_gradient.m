% Compute norm of relative gradient
function Rgrad = relative_gradient(val, x, grad, typf)
    gmax = 0.0;
    typxi = 1.0;
    for i = 1:length(grad)
        gmax = max(gmax, abs(grad(i) * max(x(i),typxi))/max(abs(val),typf));
    end
    Rgrad = gmax;
end