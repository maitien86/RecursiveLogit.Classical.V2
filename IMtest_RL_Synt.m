% Information test for RL
clear all global;
clear all;
importfile('./simulatedData/WSRL_NoLS.mat');
%  importfile('WSRL.mat');
globalVar;

notifyMail('set','amyeuphich@gmail.com','sntal2908');
isFixedUturn = true;
value = IMtest(Op.x);
result =  sprintf(' Infrmation matrix equality test (RL - synthetic): %f\n', value);
try
   notifyMail('send', result);
catch exection
   fprintf('\n Can not send email notification !!! \n');
end