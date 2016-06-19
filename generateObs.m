%   Generate Observatios with link size is include
%   The link size attributes are already computed for all pairs OD
%   --------
%   filename:   Name of file which stores the observations
%   x0:         Given parameters
%   ODpairs:    Matrix with all OD pairs
%   nbobsOD:    Number of generated obs each OD
%%
function Obs = generateObs(filename, x0, ODpairs, nbobsOD)

    global incidenceFull; 
    global Mfull;
    global Ufull;
    global Atts;
    global Op;
    global LSatt;
    global LinkSize;
    global isLinkSizeInclusive;
    
    % Generate Obs
    % ----------------------------------------------------  
    mu = 1; % MU IS NORMALIZED TO ONE
    % Parameter for the radom term
    location = 0;
    scale = 1;
    euler = 0.577215665;
    
    [lastIndexNetworkState] = size(incidenceFull, 1);
    [p q] = size(incidenceFull);
    nbobs = size(ODpairs,1);
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
    IregularNetwork = incidenceFull(1:lastIndexNetworkState,1:lastIndexNetworkState);    
    dummy = lastIndexNetworkState + 1;
    % Obs = sparse(zeros(nbobs * nbobsOD, dummy));
    % Loop over all OD pairs
    rng('shuffle');
    for n = 1:nbobs
        %n
        dest = ODpairs(n, 1);
        orig = ODpairs(n, 2);
        % get Link Size attributes
        if isLinkSizeInclusive == true            
            LinkSize = LSatt(n).value;      
            AttLc(Op.n) =  Matrix2D(LinkSize(1:lastIndexNetworkState,1:lastIndexNetworkState));
            AttLc(Op.n).Value(:,lastIndexNetworkState+1) = sparse(zeros(lastIndexNetworkState,1));
            AttLc(Op.n).Value(lastIndexNetworkState+1,:) = sparse(zeros(1, lastIndexNetworkState + 1));       
        end
        if true   
            % Get M and U
            Mfull = getM(x0, isLinkSizeInclusive); % matrix with exp utility for given beta
            M = Mfull(1:lastIndexNetworkState,1:lastIndexNetworkState);            
            addColumn = Mfull(:,dest);
            M(:,lastIndexNetworkState+1) = addColumn;
            M(lastIndexNetworkState+1,:) = zeros(1,lastIndexNetworkState+1);
            [Z, expVokBool] = getExpV(M); % vector with value functions for given beta                                                                     
            if (expVokBool == 0)
                disp('The parameters are not fesible')
                ok = false;
                return; 
            end    
            V = log(Z);
        end                
        % Get Utility
        Ufull = getU(x0, isLinkSizeInclusive); % matrix of utility for given beta
        U = Ufull(1:lastIndexNetworkState,1:lastIndexNetworkState);            
        addColumn = Ufull(:,dest);
        U(:,lastIndexNetworkState+1) = addColumn;
        U(lastIndexNetworkState+1,:) = zeros(1,lastIndexNetworkState+1);         
        addColumn = incidenceFull(:,dest);
        Incidence = IregularNetwork;
        Incidence(:,lastIndexNetworkState+1) = addColumn;
          
        % Now we have all utilities, time to generate the observations      
        for i = 1: nbobsOD
            Obs((n-1)*nbobsOD + i, 1) = dest;
            Obs((n-1)*nbobsOD + i, 2) = orig;
            k = orig;
            t = 3;
            while k ~= dummy
                ind = find(Incidence(k,:));
                nbInd = size(ind,2);
                bestUtilities = -1e6;
                for j = 1: nbInd
                    %utility = U(k,ind(j)) + V(ind(j)) + random('ev',location,scale) - euler ;                   
                    utility = U(k,ind(j)) + V(ind(j)) +  evrnd_local(location,scale) - euler ;                   
                    if bestUtilities < utility
                        bestInd = ind(j);
                        bestUtilities = utility;
                    end
                end
                if bestInd ~= dummy
                    Obs((n-1)*nbobsOD + i, t) = bestInd;
                    t = t + 1;
                end
                k = bestInd;
            end
            Obs((n-1)*nbobsOD + i, t) = dest;
        end
    end
    % Write to file::
    if isempty(filename)
        Obs = sparse(Obs);
    else
        [i,j,val] = find(Obs);
        data_dump = [i,j,val];
        save(filename,'data_dump','-ascii');
    end
    %ok = true;
    
end