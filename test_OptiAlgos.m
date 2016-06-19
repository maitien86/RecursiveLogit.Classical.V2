function [] = test_OptiAlgos(nVariable, Hessian, fixedUturn, fileObs) 

    globalVar;
    global isLinkSizeInclusive;
    global file_observations;
    global isFixedUturn;
    global RealOD;
    RealOD = spconvert(load('./Input/observationsForEstimBAI.txt'));
    RealOD = RealOD(:,1:2);
    
    tic;
    file_linkIncidence = './Input/linkIncidence.txt';
    file_AttEstimatedtime = './Input/ATTRIBUTEestimatedtime.txt';
    file_turnAngles = './Input/ATTRIBUTEturnangles.txt';
    
    if str2num(nVariable) == 4
        isLinkSizeInclusive = false;
    else
        isLinkSizeInclusive = true;
    end
    if strcmp(fixedUturn,'true') == true
        isFixedUturn = true;
    else
        isFixedUturn = false;
    end
    
    file_observations = fileObs;
    loadData;
    %nbobs = 10;
    Op = Op_structure;
    initialize_optimization_structure();
    Op.Optim_Method = OptimizeConstant.TRUST_REGION_METHOD;
    
    if strcmp(Hessian,'BHHH')    
        Op.Hessian_approx = OptimizeConstant.BHHH;
    else
        Op.Hessian_approx = OptimizeConstant.BFGS;
    end
    
    Gradient = zeros(nbobs,Op.n);
%     if isLinkSizeInclusive == false
%         Op.x = [-1.8, -1.0, -1.0]';
%     else
%         Op.x = [-1.8, -1.0, -1.0, -0.20]';
%     end

    EstimationScript;
    TXT = [nVariable,'|', Hessian,'|', fixedUturn,'|', fileObs];
    fileID = fopen('ComputationalTime_estimation.txt','at+');
    fprintf(fileID,['Estimation time|',TXT,'|Nb iterations:',num2str(Op.k),'|Elapsed time:', num2str(toc),'\n']);
    fclose(fileID);
    
    %fileID = fopen('ComputationalTime.txt','at+');
    %fprintf(fileID,['IM|',TXT,'Elapsed time:', num2str(toc),'\n']);
    %fclose(fileID);

end
