%   Get analytical Hessian matrix
%   Link size is included
%%
function [LL, grad, Hessian, Hs] = getAnalyticalHessian()
    
    global incidenceFull; 
    global Mfull;
    global Atts;
    global Obs;     % Observation
    global nbobs;  
    global LSatt;
    global LinkSize;
    global Op;
    global isLinkSizeInclusive;
    global Gradient;
    global maxstates;
    global RealOD;
  
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
    Hs = objArray(nbobs);
    gradV0 = zeros(1,Op.n);
    factor = fix(nbobs / 1832); 
    if factor < 1
        factor = 1;
    end
    % Loop over all observation
    prevOD = [0,0];
    for n = 1:nbobs
        
       %m = fix((n-1)/(factor)) + 1;
        dest = Obs(n, 1);       
        orig = Obs(n, 2);
        % get Link Size attributes
        if isLinkSizeInclusive == true
            m = find(ismember(RealOD,[dest,orig],'rows'),1);
            LinkSize = LSatt(m).value;      
            AttLc(Op.n) =  Matrix2D(LinkSize(1:lastIndexNetworkState,1:lastIndexNetworkState));
            AttLc(Op.n).Value(:,lastIndexNetworkState+1) = sparse(zeros(lastIndexNetworkState,1));
            AttLc(Op.n).Value(lastIndexNetworkState+1,:) = sparse(zeros(1, lastIndexNetworkState + 1));
        end
        if  norm([dest,orig] - prevOD) ~= 0
            %m
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
            prevOD = [dest,orig];
        end
        
        % Compute Gradient:
        %--------------------------------------------------------------
        % Get Utility
        Ufull = getU(Op.x, isLinkSizeInclusive); % matrix of utility for given beta
        U = Ufull(1:lastIndexNetworkState,1:lastIndexNetworkState);            
        addColumn = Ufull(:,dest);
        U(:,lastIndexNetworkState+1) = addColumn;
        U(lastIndexNetworkState+1,:) = zeros(1,lastIndexNetworkState+1);       
        sumInstU = 0;
        sumInstX = zeros(1,Op.n);
        seq = 2;
        a = Obs(n,seq+1); % action state after origin
        lnPn = - 1 * ((1/mu) * log(expV0));
        Gradient(n,:) = - gradV0;
        while ((a ~= 0) && (seq < (maxstates - 1)))
            k = Obs(n,seq); % current state
            if (a > (lastIndexNetworkState+1))
                a = lastIndexNetworkState + 1; % The dest index for computation is always the same
            end        
            sumInstU = sumInstU + U(k,a);          
            for i = 1 : Op.n
                sumInstX(i) = sumInstX(i) + AttLc(i).Value(k,a);
            end                 
            seq = seq + 1;
            a = Obs(n,seq+1);
        end
        Gradient(n,:) = - gradV0 + sumInstX;
        
        %-----------------------------------------------------------        
        % Pre-compute the gradient and the Hessian of M
        if mod(n-1,factor) == 0
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
        end
        Hs(n).value = H;
        Hessian = Hessian + (H - Hessian)/n;        
    end
    LL = 0;
    grad = [];
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
