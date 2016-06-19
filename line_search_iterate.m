%   A line search iterate immplemetation.
%   Tien, Mai Anh -  Diro.
%%
function Rvl = line_search_iterate()
    global Op;
    % Compute the search direction Hp = -grad
    p = linsolve(Op.H, -Op.grad);
    if p' * Op.grad > 0 
        p = -p;
    end
    % construct anonymous function for arc
    arc = @(stp) line_arc(stp,p);
    stp = OptimizeConstant.INITIAL_STEP_LENGTH;
    x = Op.x;
    grad = Op.grad;
    val = Op.value;
    fcn = @LL;
     % Get line search step 
    [x val grad stp info nfev] = line_search_asrch(fcn, x, val, grad, arc, stp , ...
                                   OptimizeConstant.NEGATIVE_CURVATURE_PARAMETER , OptimizeConstant.SUFFICIENT_DECREASE_PARAMETER, ... 
                                   OptimizeConstant.CURVATURE_CONDITION_PARAMETER ,OptimizeConstant.X_TOLERENT, ... 
                                   OptimizeConstant.MINIMUM_STEP_LENGTH, OptimizeConstant.MAXIMUM_STEP_LENGTH ,OptimizeConstant.MAX_FEV, []);
 
     info
     %	Line search successful                           
     if val <= Op.value
        Op.step = p * stp;
        Op.deltaGrad = grad - Op.grad;
        Op.deltaValue = val - Op.value;
        Op.value = val;
        Op.x = x;
        Op.grad = grad;
        % Update Hessian approxiamtion
        update_hessian_approx();
        Rvl = true;
     else
        Rvl = false;
     end
end