%
%   Trust region iterate implementation
%   MAI ANH TIEN - 9.August.2013
%
%%
function isSuccess = btr_interate()    
    global Op;
    Op.step = btr_step_steihaug_toint(Op);
    xold = Op.x;
    [newValue gradplus] = LL(Op.x + Op.step);
    Op.deltaGrad = gradplus - Op.grad;    
    % compute quadratic model
    md = Op.value + Op.grad'*Op.step + 0.5*Op.step' * Op.H * Op.step;
    rho = (Op.value - newValue) / (Op.value - md);
    if btr_accept_candidate(rho)
        Op.deltaValue = newValue - Op.value;
        Op.value = newValue;
        Op.grad = gradplus;
        update_hessian_approx();
        isSuccess = true;
    else
        isSuccess = false;
        Op.x = xold;
    end    
    % Update radius
    btr_update_radius(rho);  
end
