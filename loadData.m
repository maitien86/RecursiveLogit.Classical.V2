%   Load route choice data
%%
disp('Loading data ...')

global file_AttEstimatedtime;
global file_linkIncidence;
global file_observations;
global file_turnAngles;

%   Load link indidence
global incidenceFull;
incidenceFull = spconvert(load(file_linkIncidence));
[lastIndexNetworkState, n] = size(incidenceFull);

%   Load estimated time
%   Adjust the size to account for dummy links
%   att1 ~ EstimatedTime
global EstimatedTime;
EstimatedTime = spconvert(load(file_AttEstimatedtime));
EstimatedTime(lastIndexNetworkState,n)=0;

%   Load Left Turn
%   att2 ~ Left Turn
global TurnAngles;
TurnAngles = spconvert(load(file_turnAngles));
TurnAngles(lastIndexNetworkState,n)=0;

%   Defining u-turns (angle larger than 177 degrees)
%   att4 ~ Uturn
global Uturn;
Uturn = TurnAngles;
I = find(TurnAngles);
[nbnonzero,c] = size(I);
for i = 1:nbnonzero
    if ((Uturn(I(i)) > 3.1) || (Uturn(I(i)) < -3.1))
        Uturn(I(i)) = 1;
    else
        Uturn(I(i)) = 0;
    end
end
%   MAKE A LEFT TURN DUMMY AND SET ALL OTHER TURN INFO TO ZERO
%   leftlimit = -0.6989 % left if radians less than this (40 degrees)
%   with this critera, approx 1/5 of possible turns are considered to be left
%   LeftTurn ~ LeftTurn
global LeftTurn;
LeftTurn = TurnAngles; % create a copy of TurnAngles including all turns
for i=1:nbnonzero
    if ((TurnAngles(I(i)) < 3.1) && (TurnAngles(I(i)) > -3.1)) % NOT A U-TURN
        if (TurnAngles(I(i)) > 0) % RIGHT TURN
            TurnAngles(I(i)) = 0;
        else % Left turn
            LeftTurn(I(i)) = 0;
            % if (TurnAngles(I(i)) < -0.6989) % CREATE DUMMY, if continuous remove this if completely.
            if (TurnAngles(I(i)) < -0.5236) % CREATE DUMMY, if continuous remove this if completely.
                TurnAngles(I(i)) = -1;            
            else
                TurnAngles(I(i)) = 0;
            end
        end        
    else
        TurnAngles(I(i)) = 0; % U-TURN IS ZERO (described by att4)
        LeftTurn(I(i)) = 0;
    end
end
TurnAngles = -1 * TurnAngles; % All angles should be positive
%   Take away u-turns from LeftTurn 
I = find(LeftTurn);
[nbnonzero,c] = size(I);
for i=1:nbnonzero      
    if ((LeftTurn(I(i)) > 3.1) || (LeftTurn(I(i)) < -3.1)) % NOT A U-TURN
        LeftTurn(I(i)) = 0;
    end
end
LeftTurn = EstimatedTime;
I = find(EstimatedTime);
[nbnonzero,c] = size(I);
for i = 1:nbnonzero
    if (EstimatedTime(I(i)) > 0)
        LeftTurn(I(i)) = 1;
    else
        LeftTurn(I(i)) = 0;
    end
        
end
% For link size attribute 

global isLinkSizeInclusive;
%   Define the b vector
%   Always one at the same element
%   Always the size of the states in the "without-dummies-network" plus 1
global b;
b = zeros(lastIndexNetworkState+1,1);
b(lastIndexNetworkState+1)=1;
b_lineq = sparse(b);

%   Load observations
global Obs;
global nbobs;       % Number of observation
global maxstates;   % Number of states
global Atts;

Obs = spconvert(load(file_observations));
[nbobs, maxstates] = size(Obs);
nbobs = 10;


getAtt(); 
%% For link size attribute
% if isLinkSizeInclusive == true
%     disp('Computing Link Size attribute ... ');
%     LSatt = objArray(nbobs);
%     ok = getLinkSizeAtt();
%     if ok == false
%         fprintf(' Beta is wrong, try other parameters \n');
%     end
% end
