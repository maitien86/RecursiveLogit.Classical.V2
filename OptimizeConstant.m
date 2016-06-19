% Constants for optimization methods
classdef OptimizeConstant
    properties (Constant = true)
        NUM_ERROR = -1e-3;
        RESIDUAL  =  1e-3;
        LL_ERROR_VALUE = 99999;
        DBL_EPSILON = 2.22e-16;
        R_SR1  = 1e-8;
        
        %For line search parameters
        LINE_SEARCH_METHOD = 'LINE SEARCH METHOD';
        INITIAL_STEP_LENGTH = 1.0;
        NEGATIVE_CURVATURE_PARAMETER = 0.0;
        SUFFICIENT_DECREASE_PARAMETER = 0.0001 ;
        CURVATURE_CONDITION_PARAMETER = 0.9;
        X_TOLERENT = 2.2e-16;
        MINIMUM_STEP_LENGTH = 0;
        MAXIMUM_STEP_LENGTH = 1000;
        MAX_FEV = 10; % Maximum number of function evaluations
       
        % For trust region
        TRUST_REGION_METHOD = 'TRUST REGION METHOD';
        INITIAL_TRUST_REGION_RADIUS = 1.0;
        
        % Hessian approximation
        BFGS = 0;
        SR1 = 1;
        BHHH = 2;
        CB_BFGS = 3;
        CB_SR1 = 4;
        BIGGS = 5;    
        SSA_BFGS = 6;
        SSA_SR1 = 7;        
    end
    
    methods (Static = true)
        function name  = getHessianApprox(ind)
             switch ind
                 case 0
                     name = 'BFGS';
                 case 1
                     name = 'SR1';
                 case 2
                     name = 'BHHH';
                 case 3
                     name = 'CB-BFGS';
                 case 4
                     name = 'CB-SR1';
                 case 5 
                     name = 'BIGGS';
                 case 6
                     name = 'SSA-BFGS';
                 case 7
                     name = 'SSA-SR1';
             end
        end
    end
end