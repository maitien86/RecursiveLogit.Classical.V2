% Information test for RL
clear all global;
clear all;
importfile('./WSRL.mat');
%  importfile('WSRL.mat');
globalVar;

notifyMail('set','amyeuphich@gmail.com','sntal2908');
isFixedUturn = false;
result = IM_fulltest(Op.x);
try
   notifyMail('send', result);
catch exection
   fprintf('\n Can not send email notification !!! \n');
end