%% Compute Prediction

function results = prediction(isLS, isObs)

    globalVar;
    global TXT; 
    global SampleObs;
    notifyMail('set','amyeuphich@gmail.com','sntal2908');

    %% Data 

    file_linkIncidence = './Input/linkIncidence.txt';
    file_AttEstimatedtime = './Input/ATTRIBUTEestimatedtime.txt';
    file_turnAngles = './Input/ATTRIBUTEturnangles.txt';
    file_observations = './Input/observationsForEstimBAI.txt';
    %file_observations = './simulatedData/ObservationsAll_NoLS.txt';

    isLinkSizeInclusive = isLS;
    isFixedUturn = false;
    loadData;
    Op = Op_structure;
    initialize_optimization_structure();
    
    if isObs ==  true
        PredSample = spconvert(load('../PredSampleObs.21.40.txt'));
        note = 'OBS';
    else
        PredSample = spconvert(load('../PredSampleODs.21.40.txt'));
        note = 'ODS';
    end
    TXT = ['RL prediction:',note,':'];
    TXT = [TXT sprintf('\n Link size = %d \n', isLinkSizeInclusive)];
    nTest = round(size(PredSample,1) / 2);
    %% Estimation for tranning samples
    for i = 1: nTest
        train = PredSample(i*2-1,:);
        test  = PredSample(i*2,:);
        SampleObs = train; 
        RLestimator(train, isLS);
        LL = getPLL(test);
        TXT = [TXT sprintf('%d : %f \n', i, LL)];        
    end    
    try
            notifyMail('send', TXT);
        catch exection
            fprintf('\n Can not send email notification !!! \n');
        end
end