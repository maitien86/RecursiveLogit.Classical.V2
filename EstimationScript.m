    
    %% Optimizaing
    disp('Start Optimizing ....')
    [Op.value, Op.grad ] = LL(Op.x);
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
