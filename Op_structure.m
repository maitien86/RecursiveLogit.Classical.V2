% Optimization structure
classdef Op_structure
   properties
      Optim_Method = OptimizeConstant.LINE_SEARCH_METHOD; % or 'BTR'
      CombinedType = 'NONE';
      Hessian_approx = OptimizeConstant.BFGS;     
      n = 4;          % size of vector
      alpha = 0.0001; % Paramters for linesearch algorithm
      gamma = 0.9 ;   % Paramters for linesearch algorithm
      maxIter = 200;
      nFev = 0;       % Number of function evaluation
      tol = 0.00001;
      x = [];
      grad = [];
      deltaGrad = [];
      step = [];
      stepLength = 0;
      k = 0;
      value = 0;
      deltaValue = 0;
      H = [];   
      Ak = [];% Hessian approximation.
      Hi;
      radius = 1;   
      % trust region radius.   
      delta = 1;
      ETA1 = 0.2;
      ETA2 = 0.75;
   end
end