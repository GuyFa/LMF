x_initial=GIF.ScaffoldUVonCones;
% numDOF=(length(x)+2)/10;%82
% numDOF1=4*numDOF;%328
% numDOFScaffold=5*numDOF;%410
% numDOF1=numDOFScaffold-numDOF;
totVars=length(x_initial);%818
almostAllVars=totVars-numDOF;%736
x_initial=[x_initial([1:(numDOF1-1)], :); flip(x_initial([numDOF1:(numDOFScaffold-1)], :)); x_initial([numDOFScaffold:almostAllVars], :); flip(x_initial([(almostAllVars+1):totVars], :))];
% boundingBoxDimensionsPerRemeshing = zeros(100,9);
% boundingBoxDimensionsPerRmeshing1=1.1*boundingBoxDimensionsPerRmeshing;
iter = 0;
timesRemeshed = 0;
stepTol1=1e-3;
gradTol1=1e-3;
Etol1=1e-3;
EitersCounter=0;
numConvIters=5;
n=(length(x_initial)+length(GIF.FixedIndices))/2;
solution=zeros(length(x_initial)+length(GIF.FixedIndices),1); 
solution(setdiff(1:end,GIF.FixedIndices))=x_initial; % set free values 
solution(GIF.FixedIndices)=GIF.FixedValues; % set fixed values
solution=[solution(1:n),solution(n+1:end)];  
solution([1:numDOF1], :)=solution([1:numDOF1], :);
a=solution([(numDOF1+1):numDOFScaffold], :);
GIF.UVonCones=a;
GIF.UVonCones1=[solution([1:numDOF1], :); circshift(solution([(numDOF1+1):numDOFScaffold], :),1,1)];
GIF.UVonCones2=[solution([1:numDOF1], :); flip(solution([(numDOF1+1):numDOFScaffold], :)); a];

