% test Z
%%
   globalVar;
    % Create scale  
    %global scale;
   
    isLinkSizeInclusive = false;
    isFixedUturn = false;
    Op = Op_structure;
    initialize_optimization_structure();
    %Op.x = [-1.8 -0.933 -0.411 -3.0]';
    Op.Optim_Method = OptimizeConstant.TRUST_REGION_METHOD;
    Op.Hessian_approx = OptimizeConstant.SSA_BFGS;
    Gradient = zeros(nbobs,Op.n);
    Op.x = [-1.5,-1.5,-1.5,-1.5]';
    Mfull = getM(Op.x, false);
    [p q] = size(Mfull);
    MregularNetwork = Mfull(1:lastIndexNetworkState,1:lastIndexNetworkState);
    Ufull = getU(Op.x, false);
    UregularNetwork = Ufull(1:lastIndexNetworkState,1:lastIndexNetworkState);
    % Set LL value
    LL = 0;
    grad = zeros(1, Op.n);
    % Initialize
    expV0 = zeros(1, q + 1);
    M = MregularNetwork;
    M(:,lastIndexNetworkState+1) = sparse(zeros(lastIndexNetworkState,1));
    M(lastIndexNetworkState+1,:) = sparse(zeros(1, lastIndexNetworkState + 1));
    U = UregularNetwork;
    U(:,lastIndexNetworkState+1) = sparse(zeros(lastIndexNetworkState,1));
    U(lastIndexNetworkState+1,:) = sparse(zeros(1, lastIndexNetworkState + 1));
    gradV0 = zeros(q, Op.n);
    for i = 1:Op.n
        AttLc(i) =  Matrix2D(Atts(i).Value(1:lastIndexNetworkState,1:lastIndexNetworkState));
        AttLc(i).Value(:,lastIndexNetworkState+1) = sparse(zeros(lastIndexNetworkState,1));
        AttLc(i).Value(lastIndexNetworkState+1,:) = sparse(zeros(1, lastIndexNetworkState + 1));
    end
    dest = Obs(1, 1);
    M(1:lastIndexNetworkState ,lastIndexNetworkState + 1) = Mfull(:,dest);% addColumn;
    [n m]= size(M);
    b = zeros(n,1);
    b(n)=1;
    b_sp = sparse(b);    
    Z = ones(n,1);  
    %scale
     a = rand(n,1);
     k = 1./a;
     MI = M; MI(find(M)) = 1;
     scale = sparse((a * k') .* MI);
     
%     

    while(1)
        Z0 = sparse(exp(repmat(log(Z)',n,1) .* scale)); 
        Z1 = diag(M * Z0') + b_sp;
        if norm(log(Z) - log(Z1))<0.00001
            break;
        end
        norm(log(Z) - log(Z1))
        Z = Z1;
     
    end
    Z;