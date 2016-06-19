%   Hessian approximation
%   BHHH - BFGS - SR1 - BIGGS - Compbined_BFGS - Conbined_SR1
%%
function ok = update_hessian_approx()
    global Op;
    switch Op.Hessian_approx
        case OptimizeConstant.BHHH
            Op.H = BHHH();
            ok = true;
        case OptimizeConstant.SR1
            [Op.H ok] = SR1(Op.step, Op.deltaGrad, Op.H);
        case OptimizeConstant.BFGS
            [Op.H ok] = BFGS(Op.step, Op.deltaGrad, Op.H);            
        case OptimizeConstant.CB_BFGS
            [Op.H ok] = Combined_BFGS(Op.step, Op.deltaGrad, Op.Ak);
        case OptimizeConstant.CB_SR1
            [Op.H ok] = Combined_SR1(Op.step, Op.deltaGrad, Op.Ak);
        case OptimizeConstant.BIGGS
            [Op.H ok] = BIGGS(Op.deltaValue, Op.step, Op.deltaGrad, Op.grad, Op.H);
        case OptimizeConstant.SSA_BFGS
            [Op.H ok] = SecantStatistical_approx(Op.step, Op.deltaGrad, OptimizeConstant.SSA_BFGS);
        case OptimizeConstant.SSA_SR1
            [Op.H ok] = SecantStatistical_approx(Op.step, Op.deltaGrad, OptimizeConstant.SSA_SR1);
    end
end