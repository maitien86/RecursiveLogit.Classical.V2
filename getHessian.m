%   Get analytical Hessian matrix
%   Link size is included
%%
function Hessian = getHessian()
    
    global incidenceFull; 
    global Mfull;
    global Atts;
    global Obs;     % Observation
    global nbobs;  
    global LSatt;
    global LinkSize;
    global Op;
    global isLinkSizeInclusive;
  
    mu = 1; % MU IS NORMALIZED TO ONE
    % TO DO: compute this once and send these as parameters to this function
    [lastIndexNetworkState, nsize] = size(incidenceFull);
    [p q] = size(incidenceFull);
    expV0 = zeros(1, q + 1);
    % For the OD independence attributes
    if isLinkSizeInclusive == true
        sizeOfParams = Op.n - 1 
    else
        sizeOfParams = Op.n;
    end
      
    for i = 1 : sizeOfParams
        AttLc(i) =  Matrix2D(Atts(i).Value(1:lastIndexNetworkState,1:lastIndexNetworkState));
        AttLc(i).Value(:,lastIndexNetworkState+1) = sparse(zeros(lastIndexNetworkState,1));
        AttLc(i).Value(lastIndexNetworkState+1,:) = sparse(zeros(1, lastIndexNetworkState + 1));
    end
    
    % Initialize gradient and hessian of matrix M
    gradM = objArray(Op.n);
    hessianM = objMatrix(Op.n, Op.n);
    hessianZ0 = zeros(Op.n);
    Hessian = zeros(Op.n);
    H = zeros(Op.n);
    gradV0 = zeros(1,Op.n);
    % Loop over all observation
    for n = 1:nbobs
        n
        dest = Obs(n, 1);       
        orig = Obs(n, 2);
        % get Link Size attributes
        if isLinkSizeInclusive == true
            LinkSize = LSatt(n).value;      
            AttLc(Op.n) =  Matrix2D(LinkSize(1:lastIndexNetworkState,1:lastIndexNetworkState));
            AttLc(Op.n).Value(:,lastIndexNetworkState+1) = sparse(zeros(lastIndexNetworkState,1));
            AttLc(Op.n).Value(lastIndexNetworkState+1,:) = sparse(zeros(1, lastIndexNetworkState + 1));
        end
        if true %(expV0(dest) == 0)   
            % Get M and U
            Mfull = getM(Op.x, isLinkSizeInclusive); % matrix with exp utility for given beta
            M = Mfull(1:lastIndexNetworkState,1:lastIndexNetworkState);            
            addColumn = Mfull(:,dest);
            M(:,lastIndexNetworkState+1) = addColumn;
            M(lastIndexNetworkState+1,:) = zeros(1,lastIndexNetworkState+1);
            M = sparse(M);
            [expV, expVokBool] = getExpV(M); % vector with value functions for given beta                                                                     
            if (expVokBool == 0)
                LL = OptimizeConstant.LL_ERROR_VALUE;
                grad = [1 1 1 1 1]';
                disp('The parameters are not fesible')
                return; 
            end  
            % expV ~ Z           
            expV0(dest) = expV(Obs(n,2)); % Get the exp(V0)
            gradExpV = getGradExpV(M, AttLc, Op, expV);
            for i = 1:Op.n
                gradV0(i) =  gradExpV(i).Value(orig)/expV(orig);
            end 
        end
        % Pre-compute the gradient and the Hessian of M
        I = speye(size(M));  
        A = I - M;
        for i = 1: Op.n
            gradM(i).value =  M .* (AttLc(i).Value);
        end
        for i = 1: Op.n
            for j = 1: Op.n
                u = AttLc(i).Value .* AttLc(j).Value;
                v = M .* u;
                hessianM(i,j).value = v;
            end
        end
        % Compute the hessian of expV
        for i = 1: Op.n
            for j = 1: Op.n
                u = hessianM(i,j).value * expV;
                p1 = (A\u);
                u = gradM(i).value * gradExpV(j).Value;
                p2 = (A\u);
                u = gradM(i).value * expV;
                v = gradM(j).value * (A\u);
                p3 = A\v;
                p = (p1 + p2 + p3);
                hessianZ0(i,j) = p(orig);
            end
        end
        % Compute Hessian of LL
        for i = 1: Op.n
            for j = 1: Op.n
                H(i,j) = (hessianZ0(i,j)/expV(orig) - gradV0(i) * gradV0(j));
            end
        end
        Hessian = Hessian + (H - Hessian)/n;        
    end
end

function [gradExpV] = getGradExpV(M, Att, op, expV)
    I = speye(size(M));  
    A = I - M; 
    for i = 1:op.n
        u = M .* (Att(i).Value); 
        v = u * expV; 
        gradExpV(i) = Matrix2D(A\v); 
    end
end
