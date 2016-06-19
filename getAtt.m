%   Get attribute 
function [] = getAtt()

    global incidenceFull;
    global EstimatedTime;
    global TurnAngles;
    global LeftTurn;
    global Uturn;
    global Op;
    global Atts;
    global nbobs;
    global attsType;
    global isLinkSizeInclusive;
    global LSatt;
    Atts  = objArray(Op.n);
    attsType = [0,0,0,0,1]; % 0 if the correspnding attribute is observation independent, 1 otherwise   
    Incidence = incidenceFull;
    %% Observation independent attributes
    Atts(1).value = (Incidence .* EstimatedTime);
    Atts(2).value = (Incidence .* TurnAngles);
    Atts(3).value = (Incidence .* LeftTurn);
    Atts(4).value = (Incidence .* Uturn);
    %% Observation specific attributes - Link Size
    
    Atts(5).value = objArray(nbobs);
    if isLinkSizeInclusive == true
      disp('Computing Link Size attribute ... ');
      LSatt = objArray(nbobs);
      ok = getLinkSizeAtt();
      if ok == false
         fprintf(' Beta is wrong, try other parameters \n');
      end
    end
    for i = 1:nbobs
        Atts(5).value = LSatt;
    end
end