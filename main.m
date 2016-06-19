....
%   Linked based network route choice model with unrestricted choice set
%   Optmization algorithm
%   MAI ANH TIEN - DIRO
%   29 - July - 2013
%   MAIN PROGRAM
%   ---------------------------------------------------
%%
Credits;
% Initialize email notification
global resultsTXT; 
notifyMail('set','amyeuphich@gmail.com','sntal2908');

clear all global;
clear all;
%importfile('./simulatedData/WSPSL.mat');
importfile('WSRL.mat');
globalVar;

%------------------------------------------------------
%   Load data
global file_AttEstimatedtime;
global file_linkIncidence;
global file_observations;
global file_turnAngles;
global nbobs;       % Number of observation
global Op;
global isLinkSizeInclusive;
global Gradient;

isLinkSizeInclusive = false;
file_linkIncidence = './Input/linkIncidence.txt';
file_AttEstimatedtime = './Input/ATTRIBUTEestimatedtime.txt';
file_turnAngles = './Input/ATTRIBUTEturnangles.txt';
file_observations = './Input/observationsForEstimBAI.txt';
%file_observations = './simulatedData/Observations1000.txt';

%loadData;

% Initialize the optimizer structure
Op = Op_structure;
initialize_optimization_structure();
Op.x = [-2.494 -0.933 -0.411 -4.459]';
Op.Optim_Method = OptimizeConstant.TRUST_REGION_METHOD;
Op.Hessian_approx = OptimizeConstant.BFGS;
%Op.Optim_Method = OptimizeConstant.LINE_SEARCH_METHOD;
Gradient = zeros(nbobs,Op.n);

% Generate Observations
% createSimulatedObs;

%---------------------------
%Starting optimization
tic ;
%progTest
disp('Start Optimizing ....')
[Op.value, Op.grad ] = getLL();
PrintOut(Op);
% print result to string text
header = [sprintf('%s \n',file_observations) Op.Optim_Method];
header = [header sprintf('\nNumber of observations = %d \n', nbobs)];
header = [header sprintf('Hessian approx methods = %s \n', OptimizeConstant.getHessianApprox(Op.Hessian_approx))];
resultsTXT = header;
%------------------------------------------------
while (true)    
  Op.k = Op.k + 1;
  if strcmp(Op.Optim_Method,OptimizeConstant.LINE_SEARCH_METHOD);
    ok = line_search_iterate();
    if ok == true
        PrintOut(Op);
    else
        disp(' Unsuccessful process ...')
        break;
    end
  else
    ok = btr_interate();
    PrintOut(Op);
  end
  [isStop, Stoppingtype, isSuccess] = CheckStopping(Op);  
  %----------------------------------------
  % Check stopping criteria
  if(isStop == true)
      isSuccess
      fprintf('The algorithm stops, due to %s', Stoppingtype);
      resultsTXT = [resultsTXT sprintf('The algorithm stops, due to %s \n', Stoppingtype)];
      break;
  end
end
%   Compute variance - Covariance matrix
%PrintOut(Op);
getCov;
%   Finishing ...
ElapsedTtime = toc
resultsTXT = [resultsTXT sprintf('\n Number of function evaluation %d \n', Op.nFev)];
resultsTXT = [resultsTXT sprintf('\n Estimated time %d \n', ElapsedTtime)];
try
   notifyMail('send', resultsTXT);
catch exection
   fprintf('\n Can not send email notification !!! \n');
end

