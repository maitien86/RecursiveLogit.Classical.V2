%   Calculate variance - covariance matrix
%%
global nbobs;
global resultsTXT;
global Op;
fprintf('\n Getting analytical Hessian ....\n');
Hessian = getHessian();
InvH = inv(Hessian);
CoV = (InvH * BHHH() * InvH)/nbobs;
resultsTXT = [resultsTXT sprintf('\n Hessian is: \n') sprintf([repmat('%f\t', 1, size(Hessian, 2)) '\n'], Hessian')];
resultsTXT = [resultsTXT sprintf('\n Variance - Covariance is: \n') sprintf([repmat('%f\t', 1, size(CoV, 2)) '\n'], CoV')];
StandarError = zeros(1,Op.n);
for i = 1:Op.n
  StandarError(i) = sqrt(CoV(i,i));
end
resultsTXT = [resultsTXT sprintf('\n Standard error is: \n') sprintf(' %f ',StandarError)];