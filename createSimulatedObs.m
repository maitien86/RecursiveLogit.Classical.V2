%   Generate obs
%%
global Op;
global Obs;
globalVar;

disp('Observation generating ....')
Op.n = 4;
isLinkSizeInclusive = false;
isFixedUturn = false;
x0 = [-2.0, -1.0, -1.0, -20.0, -0.20 ]';
% ODpairs = Obs(:,1:2);
% For tutorial
% ODpairs = [27,1;27,2];
% x0 = [-2.0;0;0;0];
%------------------------------
nbobsOD = 1;
filename = './ExampleTutorial/SyntheticObs.txt';
generateObs(filename, x0, ODpairs, nbobsOD);
