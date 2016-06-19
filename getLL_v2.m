% Compute the loglikelohood value and its gradient.
%-----------------------------------------------------
function [LL, grad] = getLL_v2()

    global incidenceFull; 
    global Gradient;
    global Op;
    global Mfull;
    global Ufull;
    global Atts;
    global Obs;     % Observation
    global nbobs;  
    global SampleObs;
    global attsType;
    
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
    AttsLc = objArray(Op.n); 
    %------------------------
    for t = 1:nbobs
        dest = Obs(sample(t), 1);
        orig = Obs(sample(t), 2);
        for i = 1: Op.n
            if attsType(i) == 0
               AttsLc(i).value = Atts(i).value;
            else
               AttsLc(i).value = Atts(i).value(sample(t)).value;
            end
        end
        % Get M and U
        Mfull = getM(Op.x,AttsLc); % matrix with exp utility for given beta
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
        gradV0 = getGradV0(M, AttsLc, Op, Z, orig);                            
        % Get Utility
        Ufull = getU(Op.x,AttsLc); % matrix of utility for given beta
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
                sumInstX(j) = sumInstX(j) + AttsLc(j).value(path(i),mIndex);
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