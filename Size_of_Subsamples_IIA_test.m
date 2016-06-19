function [] = Size_of_Subsamples_IIA_test(nVariable, fileObs) 

    globalVar;
    global isLinkSizeInclusive;
    global file_observations;
    global isFixedUturn;
    tic;
    file_linkIncidence = './Input/linkIncidence.txt';
    file_AttEstimatedtime = './Input/ATTRIBUTEestimatedtime.txt';
    file_turnAngles = './Input/ATTRIBUTEturnangles.txt';
    % Get real OD
    global RealOD;
    RealOD = spconvert(load('./Input/observationsForEstimBAI.txt'));
    RealOD = RealOD(:,1:2);
    
    if str2num(nVariable) == 4
        isLinkSizeInclusive = false;
    else
        isLinkSizeInclusive = true;
    end
    isFixedUturn = true;
    %isFixedUturn = false;
    file_observations = fileObs;
    
    loadData;
    
   % nbobs
    %% Optimizing with sub choice set
    incidenceFull = spconvert(load('./Input/linkIncidence_sub1000.txt'));
    LSatt;
    % Resample Obs
    index = ones(size(Obs,1),1);
    for i = 1:size(Obs,1)
        path = Obs(i,:);
        if isPathExiste(path) == false
            index(i) = 0;
        end
    end
    i = 1;
    while(i <= size(Obs,1))
        path = Obs(i,:);
        if isPathExiste(path) == false
            Obs(i,:) = [];
            i = i-1;
        end
        i = i+1;
    end
%     t = 1;
%     for i = 1:size(index,1)
%         if index(i) == 1
%             LSatt(t).value = LSatt(i).value;
%             t = t+1;
%         end
%     end    
    nbobs = size(Obs,1);
    nbobs
    
end
function [ok] = isPathExiste(path)
    global incidenceFull;
    ok = true;
    for j=2: size(find(path),2)-1
        if incidenceFull(path(j), path(j+1)) == 0
            ok = false;
            return;
        end
    end

end
