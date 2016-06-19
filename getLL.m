% Compute the loglikelohood value and its gradient.
%%
function [LL, grad] = getLL()

    global incidenceFull; 
    global Gradient;
    global Op;
    global Mfull;
    global Ufull;
    global Atts;
    global Obs;     % Observation
    global nbobs;  
    global isLinkSizeInclusive;
    global SampleObs;
    
    if isempty(SampleObs)
        SampleObs = 1:nbobs;
    end
    sample = SampleObs;
    % Get the LL value - How to get it:
    % ----------------------------------------------------
    % If Link size is included
    if (isLinkSizeInclusive == true)
        [LL, grad] = getODspecLL();
        return;
    end
    mu = 1; % MU IS NORMALIZED TO ONE
    % TO DO: compute this once and send these as parameters to this function
    [lastIndexNetworkState, nsize] = size(incidenceFull);    
    Mfull = getM(Op.x, false);
    [p q] = size(Mfull);
    MregularNetwork = Mfull(1:lastIndexNetworkState,1:lastIndexNetworkState);
    Ufull = getU(Op.x, false);
    UregularNetwork = Ufull(1:lastIndexNetworkState,1:lastIndexNetworkState);
    Mfull = sparse(Mfull);
    Ufull = sparse(Ufull);
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
    U = sparse(U);
    M = sparse(M);
    gradV0 = zeros(q, Op.n);
    for i = 1:Op.n
        AttLc(i) =  Matrix2D(Atts(i).Value(1:lastIndexNetworkState,1:lastIndexNetworkState));
        AttLc(i).Value(:,lastIndexNetworkState+1) = sparse(zeros(lastIndexNetworkState,1));
        AttLc(i).Value(lastIndexNetworkState+1,:) = sparse(zeros(1, lastIndexNetworkState + 1));
    end
    % LOOp ON ALL OBSERVATIONS
    % Compute the LL and gradient.
    for t = 1:nbobs   
        dest = Obs(sample(t), 1);
        if true %(expV0(dest) == 0)   
            M(1:lastIndexNetworkState ,lastIndexNetworkState + 1) = Mfull(:,dest);% addColumn;
            [expV, expVokBool] = getExpV(M); % vector with value functions for given beta                                                                     
            if (expVokBool == 0)
                LL = OptimizeConstant.LL_ERROR_VALUE;
                grad = ones(Op.n,1);
                disp('The parameters not fesible')
                return; 
            end                      
            expV0(dest) = expV(Obs(sample(t),2)); % Get the exp(V0)
            gradV0(dest,:) = getGradV0(M, AttLc, Op, expV, Obs(sample(t),2));
        else
        end
        sumInstU = 0;
        sumInstX = zeros(1,Op.n);
        % Value function at origin
        lnPn = - 1 * ((1/mu) * log(expV0(dest)));
        Gradient(t,:) = - gradV0(dest,:);
        path = Obs(sample(t),:);
        lpath = size(find(path),2);
        for i = 2:lpath - 1
            sumInstU = sumInstU + Ufull(path(i),path(i+1)) ;
            for j = 1:Op.n
                sumInstX(j) = sumInstX(j) + Atts(j).Value(path(i),path(i+1));
            end          

        end         
        Gradient(t,:) = - gradV0(dest,:) + sumInstX;
        lnPn = lnPn + ((1/mu)*sumInstU) ;  
        LL =  LL + (lnPn - LL)/t;
        grad = grad + (Gradient(t,:) - grad)/t;
        Gradient(t,:) = - Gradient(t,:);
    end
    LL = -1 * LL; % IN ORDER TO HAVE A MIN PROBLEM
    grad =  - grad';
end

%%
% Compute loglikelihood value with link size attribute 
%-----------------------------------------------------
function [LL, grad] = getODspecLL()

    global incidenceFull; 
    global Gradient;
    global Op;
    global Mfull;
    global Ufull;
    global Atts;
    global Obs;     % Observation
    global nbobs;  
    global LSatt;
    global LinkSize;
    global SampleObs;
    global RealOD;
    
    if isempty(SampleObs)
        SampleObs = 1:nbobs;
    end
    sample = SampleObs;
    sample =1:nbobs;
    
    % Get the LL value
    % ----------------------------------------------------
  
    mu = 1; % MU IS NORMALIZED TO ONE
    % TO DO: compute this once and send these as parameters to this function
    [lastIndexNetworkState, nsize] = size(incidenceFull);
    [p q] = size(incidenceFull);
    LL = 0;
    grad = zeros(1, Op.n);
    % For the OD independence attributes
    for i = 1 : Op.n - 1
        AttLc(i) =  Matrix2D(Atts(i).Value(1:lastIndexNetworkState,1:lastIndexNetworkState));
        AttLc(i).Value(:,lastIndexNetworkState+1) = sparse(zeros(lastIndexNetworkState,1));
        AttLc(i).Value(lastIndexNetworkState+1,:) = sparse(zeros(1, lastIndexNetworkState + 1));
    end
    % LoOp over all observation
    for t = 1:nbobs
        dest = Obs(sample(t), 1);
        orig = Obs(sample(t), 2);
        % get Link Size attributes
        %getLinkSizeAtt(); % For test only
        m = find(ismember(RealOD,[dest,orig],'rows'),1);
        LinkSize = LSatt(m).value;       
        AttLc(Op.n) =  Matrix2D(LinkSize(1:lastIndexNetworkState,1:lastIndexNetworkState));
        AttLc(Op.n).Value(:,lastIndexNetworkState+1) = sparse(zeros(lastIndexNetworkState,1));
        AttLc(Op.n).Value(lastIndexNetworkState+1,:) = sparse(zeros(1, lastIndexNetworkState + 1));       
        if true   
            % Get M and U
            Mfull = getM(Op.x, true); % matrix with exp utility for given beta
            M = Mfull(1:lastIndexNetworkState,1:lastIndexNetworkState);            
            addColumn = Mfull(:,dest);
            M(:,lastIndexNetworkState+1) = addColumn;
            M(lastIndexNetworkState+1,:) = zeros(1,lastIndexNetworkState+1);
            [Z, expVokBool] = getExpV(M); % vector with value functions for given beta                                                                     
            if (expVokBool == 0)
                LL = OptimizeConstant.LL_ERROR_VALUE;
                grad = ones(Op.n,1);
                disp('The parameters are not fesible')
                return; 
            end    
            expV0 = Z(orig);
            gradV0 = getGradV0(M, AttLc, Op, Z, orig);            
        end                
        % Get Utility
        Ufull = getU(Op.x, true); % matrix of utility for given beta
        U = Ufull(1:lastIndexNetworkState,1:lastIndexNetworkState);            
        addColumn = Ufull(:,dest);
        U(:,lastIndexNetworkState+1) = addColumn;
        U(lastIndexNetworkState+1,:) = zeros(1,lastIndexNetworkState+1);       
        sumInstU = 0;
        sumInstX = zeros(1,Op.n);
        seq = 2;
        a = Obs(sample(t),seq+1); % action state after origin
        lnPn = - 1 * ((1/mu) * log(expV0));
        Gradient(t,:) = - gradV0;
        path = Obs(sample(t),:);
        lpath = size(find(path),2);
        for i = 2:lpath - 1
            mIndex = min(path(i+1), lastIndexNetworkState +1);
            sumInstU = sumInstU + U(path(i),mIndex) ;
            for j = 1:Op.n
                sumInstX(j) = sumInstX(j) + AttLc(j).Value(path(i),mIndex);
            end          

        end         
        Gradient(t,:) = - gradV0 + sumInstX;
        lnPn = lnPn + ((1/mu)*sumInstU) ;  
        LL =  LL + (lnPn - LL)/t;
        grad = grad + (Gradient(t,:) - grad)/t;
        Gradient(t,:) = - Gradient(t,:);
    end
    LL = -1 * LL; % IN ORDER TO HAVE A MIN PROBLEM
    grad =  - grad';
end