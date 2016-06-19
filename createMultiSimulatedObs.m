%   Generate obs
%%
global Op;
global Obs;
globalVar;
ODpairs = Obs(:,1:2);
disp('Observation generating ....')
Op.n = 5;
isLinkSizeInclusive = true;
isFixedUturn = false;
x0 = [-2.0, -1.0, -1.0, -20.0, -1.0 ]';
%x0 = [-2.0, -1.0, -1.0, -20.0, -0.2 ]';
%------------------------------
nbobsOD = 10;
folder  = './MCSyntheticObs18320/';
for i = 1:20
    i
    filename = ['./MCSyntheticObs18320/WithLS',num2str(i),'.txt'];%'./ExampleTutorial/SyntheticObs.txt';
    generateObs(filename, x0, ODpairs, nbobsOD);
end

%filename = ['./MCSyntheticObs1000/NoLS',num2str(0),'.txt'];%'./ExampleTutorial/SyntheticObs.txt';
%cutODpairs = ODpairs(1:10,:);
%generateObs(filename, x0, cutODpairs, nbobsOD);


Op.n = 4;
isLinkSizeInclusive = false;
isFixedUturn = false;
x0 = [-2.0, -1.0, -1.0, -20.0]';
%x0 = [-2.0, -1.0, -1.0, -20.0, -0.2 ]'; ------------------------------
nbobsOD = 1;
folder  = './MCSyntheticObs18320/';
for i = 1:20
    i
    filename = ['./MCSyntheticObs18320/NoLS',num2str(i),'.txt'];%'./ExampleTutorial/SyntheticObs.txt';
    generateObs(filename, x0, ODpairs, nbobsOD);
end