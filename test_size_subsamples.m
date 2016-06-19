folder  = './MCSyntheticObs/';
for i = 1:20
    filename = ['./MCSyntheticObs/WithLS',num2str(i),'.txt'];%'./ExampleTutorial/SyntheticObs.txt';
    Size_of_Subsamples_IIA_test('4',filename);
end