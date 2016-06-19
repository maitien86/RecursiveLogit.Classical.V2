%
%   Print iterative result to screen
%%
function [] = PrintOut(Op)
      global resultsTXT; 
      fprintf('[Iteration]: %d\n', Op.k);
      fprintf('     LL = %f\n', Op.value);
      fprintf('     x = \n');
      fprintf('         %i\n', Op.x');
      fprintf('     norm of step = %f\n', norm(Op.step));
      fprintf('     radius = %f\n', Op.delta);  
      fprintf('     Norm of grad = %f\n', norm(Op.grad));
      relatice_grad = relative_gradient(Op.value, Op.x, Op.grad, 1.0);
      fprintf('     Norm of relative gradient = %f\n', relatice_grad);
      fprintf('     Number of function evaluation = %f\n', Op.nFev);
      
      % To string
      resultsTXT = [resultsTXT sprintf('[Iteration]: %d\n', Op.k)];
      resultsTXT = [resultsTXT sprintf('     LL = %f\n', Op.value)];
      resultsTXT = [resultsTXT sprintf('     x = \n')];
      resultsTXT = [resultsTXT sprintf('         %i\n', Op.x')];
      resultsTXT = [resultsTXT sprintf('     norm of step = %f\n', norm(Op.step))];
      resultsTXT = [resultsTXT sprintf('     radius = %f\n', Op.delta)];  
      resultsTXT = [resultsTXT sprintf('     Norm of grad = %f\n', norm(Op.grad))];
      resultsTXT = [resultsTXT sprintf('     Norm of relative gradient = %f\n', relatice_grad)];    
      resultsTXT = [resultsTXT sprintf('     Number of function evaluation = %d\n', Op.nFev)]; 
end