% Compute the loglikelohood value and its gradient.
%%
function [LL, grad] = getPLL(sample)

    global incidenceFull; 
    global Gradient;
    global Op;
    global Mfull;
    global Ufull;
    global Atts;
    global Obs;     % Observation
    global nbobs;  
    global maxstates;
    global isLinkSizeInclusive;
    % Get the LL value - How to get it:
    % ----------------------------------------------------
    % If Link size is included
    if (isLinkSizeInclusive == true)
        [LL, grad] = getODspecPLL(sample);
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
    for i = 1:Op.n
        AttLc(i) =  Matrix2D(Atts(i).Value(1:lastIndexNetworkState,1:lastIndexNetworkState));
        AttLc(i).Value(:,lastIndexNetworkState+1) = sparse(zeros(lastIndexNetworkState,1));
        AttLc(i).Value(lastIndexNetworkState+1,:) = sparse(zeros(1, lastIndexNetworkState + 1));
    end
    % LOOp ON ALL OBSERVATIONS
    % Compute the LL and gradient.
    nobs = size(find(sample),2);
    for t = 1:nobs      
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
        else
        end
        sumInstU = 0;
        lnPn = - 1 * ((1/mu) * log(expV0(dest)));
        path = Obs(sample(t),:);
        lpath = size(find(path),2);
        % Compute regular attributes
        for i = 2:lpath - 1
            sumInstU = sumInstU + Ufull(path(i),path(i+1)) ;
        end
        lnPn = lnPn + ((1/mu)*sumInstU) ;  
        LL =  LL + lnPn;
    end
    LL = 1 * LL; % IN ORDER TO HAVE A MIN PROBLEM
end

%%
% Compute loglikelihood value with link size attribute 
%-----------------------------------------------------
function [LL, grad] = getODspecPLL(sample)

    global incidenceFull; 
    global Op;
    global Mfull;
    global Ufull;
    global Atts;
    global Obs;     % Observation
    global maxstates;
    global LSatt;
    global LinkSize;
    
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
    nobs = size(find(sample),2);
    for t = 1:nobs
        dest = Obs(sample(t), 1);
        orig = Obs(sample(t), 2);
        % get Link Size attributes
        %getLinkSizeAtt(); % For test only
        LinkSize = LSatt(sample(t)).value;      
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
        end                
        % Get Utility
        Ufull = getU(Op.x, true); % matrix of utility for given beta
        U = Ufull(1:lastIndexNetworkState,1:lastIndexNetworkState);            
        addColumn = Ufull(:,dest);
        U(:,lastIndexNetworkState+1) = addColumn;
        U(lastIndexNetworkState+1,:) = zeros(1,lastIndexNetworkState+1);       
        sumInstU = 0;
        lnPn = - 1 * ((1/mu) * log(expV0));
        
        path = Obs(sample(t),:);
        lpath = size(find(path),2);
        % Compute regular attributes
        for i = 2:lpath - 1
            sumInstU = sumInstU + Ufull(path(i),path(i+1)) ;
        end
        lnPn = lnPn + ((1/mu)*sumInstU) ;
        LL =  LL + lnPn ;
    end
end