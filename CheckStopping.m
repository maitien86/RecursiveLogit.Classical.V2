function [isStop, Stoppingtype, isSuccess] = CheckStopping(op)
isStop = true;
if op.k > op.maxIter
    Stoppingtype = 'NUMBER OF ITERATION';
    isSuccess = false;
    return;
end;
if norm(op.grad) < op.tol
    Stoppingtype = 'GRADIENT';
    isSuccess = true;
    return;
end;
if norm(op.step) < op.tol * op.tol
    Stoppingtype = 'TOO SMALL STEP';
    isSuccess = false;
    return;
end;
relatice_grad = relative_gradient(op.value, op.x, op.grad, 1.0);
if relatice_grad < op.tol
    Stoppingtype = 'RELATIVE GRADIENT';
    isSuccess = true;
    return;
end;
isStop = false;
Stoppingtype = 'NONE';
isSuccess = false;
end
