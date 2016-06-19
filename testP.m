% test P
%%
    global Obs;
    beta = [-1.8,-0.9,-0.8,-4.0,-0.5]';
    loadData; % Recursive logit
    Mfull = getM(beta,true);    
    M = Mfull(1:lastIndexNetworkState,1:lastIndexNetworkState); 
    M(:,lastIndexNetworkState+1) = sparse(zeros(lastIndexNetworkState,1));
    M(lastIndexNetworkState+1,:) = sparse(zeros(1, lastIndexNetworkState + 1));  
    
    for n = 1:nbobs
        n
        dest = Obs(n, 1);
        orig = Obs(n, 2);      
        % Get probabilities       
        M(1:lastIndexNetworkState ,lastIndexNetworkState + 1) = Mfull(:,dest);
        [expV, expVokBool] = getExpV(M); % vector with value functions for given beta                                                                     
        if (expVokBool == 0)
           isDone = false;
           disp('ExpV is not fesible')
           return; 
        end  
        P = getP(expV, M);
        Q = find(P > 1.5);
        if (~isempty(Q))
            disp('Out');
            break;
        end
        
    end