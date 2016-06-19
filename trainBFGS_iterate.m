function ok = trainBFGS_iterate()
    global Op;
    % Compute the search direction Hp = -grad
    p = linsolve(Op.H, -Op.grad);
    if p' * Op.grad > 0 
        p = -p;
    end
    % construct anonymous function for arc
    x = Op.x;
    step = 1.0;
    [f1 g1] = LL(x + step * p);
    if f1 < Op.value
        while true
            step = step * 2;
            f = f1;
            g = g1;
            [f1 g1] = LL(x + step * p);
            if f1 > Op.value
                step = step/2;
                break;
            end
        end
    else
        while true
            step = step/2;
            [f1 g1] = LL(x + step * p);
            if f1 <= Op.value
                break;
            end
        end
        f = f1;
        g = g1;
    end
    Op.step = p * step;
    Op.deltaGrad = g - Op.grad;
    Op.deltaValue = f - Op.value;
    Op.value = f;
    Op.x = x + p * step;
    Op.grad = g;
    % Update Hessian approxiamtion
    update_hessian_approx();
    ok = true;
end