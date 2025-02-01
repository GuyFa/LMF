if iter ~= 0
    x_initial=GIF.ScaffoldUVonCones;
    x_initial=[x_initial([1:(numDOF1-1)], :); flip(x_initial([numDOF1:(numDOFScaffold-1)], :)); x_initial([numDOFScaffold:almostAllVars], :); flip(x_initial([(almostAllVars+1):totVars], :))];

    stepTol1=max(1e-12, 1e-3 * rate);
    gradTol1=max(1e-12, 1e-3 * rate);
    Etol1=max(1e-12, 1e-3 * rate);

    EitersCounter=0;
    numConvIters=5;
end

use_GPU = false;

if size(GIF.F1, 1) > 5e5
    use_GPU = true;
end

n=(length(x_initial)+length(GIF.FixedIndices))/2;
SPDHessian=true;

temp_scaling = 1;
%we set the following arrays for the case that cholesky fails 3 times and
%avoid problems while building the scaffold
solution=zeros(length(x_initial)+length(GIF.FixedIndices),1);
solution(setdiff(1:end,GIF.FixedIndices))=x_initial; % set free values
solution(GIF.FixedIndices)=GIF.FixedValues; % set fixed values
solution=[solution(1:n),solution(n+1:end)];
solution([1:numDOF1], :)=solution([1:numDOF1], :);

a=solution([(numDOF1+1):numDOFScaffold], :);
GIF.UVonCones=a;
GIF.UVonCones1=[solution([1:numDOF1], :); circshift(solution([(numDOF1+1):numDOFScaffold], :),1,1)];
GIF.UVonCones2=[solution([1:numDOF1], :); flip(solution([(numDOF1+1):numDOFScaffold], :)); a];

% run Newton
[foldoversIndices, fz, fbz ,E]=symmetricDirichletEnergy(x_initial,GIF.TJ_fz,GIF.TJ_fbz,GIF.TArea,GIF.FixedIndices,GIF.FixedValues,SPDHessian, 0);
GIF.initTotalEnergy=E;
GIF.initEnergy=E;
timesRemeshed = timesRemeshed + 1;
EitersCounter=0;
if use_GPU
    % move arrays to gpu
    J_fz_gpu=gpuArray(GIF.TJ_fz);
    J_fbz_gpu=gpuArray(GIF.TJ_fbz);
    Area_gpu=gpuArray(GIF.TArea);
    FixedValues_gpu=gpuArray(GIF.FixedValues);
end

for i=1:100
    Eprev=E;

    % ****find energy,gradient, hessian****
    if use_GPU
        x_gpu=gpuArray(x_initial);
        [foldoversIndices, fz_gpu,fbz_gpu,E_gpu,G_gpu,H_gpu ] = symmetricDirichletEnergy(x_gpu,J_fz_gpu,J_fbz_gpu,Area_gpu,GIF.FixedIndices,FixedValues_gpu,SPDHessian, 0);

        % move to cpu
        E=gather(E_gpu);
        G=gather(G_gpu);
        H=gather(H_gpu);
        fz=gather(fz_gpu);
        fbz=gather(fbz_gpu);
    else
        [foldoversIndices, fz,fbz,E,G,H ] = symmetricDirichletEnergy(x_initial,GIF.TJ_fz,GIF.TJ_fbz,GIF.TArea,GIF.FixedIndices,GIF.FixedValues,SPDHessian, 0);
    end
    GIF.foldoversIndices1 = foldoversIndices;

    % ****solve****
    H1=[H([numDOF1:(numDOFScaffold-1)], [numDOF1:(numDOFScaffold-1)]), H([numDOF1:(numDOFScaffold-1)], [(almostAllVars+1):totVars]); H([(almostAllVars+1):totVars], [numDOF1:(numDOFScaffold-1)]), H([(almostAllVars+1):totVars], [(almostAllVars+1):totVars])];
    G1=[G([numDOF1:(numDOFScaffold-1)], :);G([(almostAllVars+1):totVars], :)];

    [R,p]=chol(H1);
    stepTol=stepTol1;
    gradTol=stepTol1;
    Etol=stepTol1;

    if p ~= 0
        [R,p]=chol(H1 + eye(size(H1)));
        stepTol=1e-3 * stepTol1;
        gradTol=1e-3 * stepTol1;
        Etol=1e-3 * stepTol1;
    end

    if p ~= 0
        [R,p]=chol(eye(size(H1)));
    end
    if p ~= 0
        GIF = rmfield(GIF,{'SF', 'SV', 'SJ_fz', 'SJ_fbz', 'SArea','TArea', 'EditedJ_fz', 'EditedJ_fbz','TJ_fz', 'TJ_fbz'});
        clearvars -except GIF x_initial numR areaOriginal SPDHessian numDOF1 totVars almostAllVars iter holeCoor numFixedDOFScaffold timesRemeshed arr1 p gradTol1 stepTol1 Etol1 EitersCounter numConvIters terminationCondition rate lambda E LAST_EVAL_TO_COUT
    end
    if p ~=0
        return;
    end
    d1=R\(R'\(-G1));

    d=zeros(totVars, 1);
    d(numDOF1:(numDOFScaffold-1))=d1(1:numDOF);
    d((almostAllVars+1):totVars)=d1((numDOF+1):(2*numDOF));

    % ****line search****
    d_full=zeros(length(d)+length(GIF.FixedIndices),1);
    d_full(setdiff(1:end,GIF.FixedIndices))=d; % set free values
    d_full(GIF.FixedIndices)=0; % set all fixed values to 0

    tJ_fz1 = GIF.TJ_fz;
    tJ_fz1(foldoversIndices, :) = [];
    tJ_fbz1 = GIF.TJ_fbz;
    tJ_fbz1(foldoversIndices, :) = [];
    tArea = GIF.TArea;
    tArea(foldoversIndices, :) = [];
    d_complex=complex(d_full(1:n),d_full(n+1:end));
    d_fz=tJ_fz1*d_complex;
    d_fbz=tJ_fbz1*d_complex;

    % locally injective line search
    tMax = lineSearchLocalInjectivity( fz,fbz,d_fz,d_fbz);

    % decreasing energy line search
    t = lineSearchDecreasingEnergy( fz,fbz,d_fz,d_fbz,tArea,E,G,d,tMax );
    % ****do step****

    x_initial=x_initial+t*d;
    iter=iter+1;

    temp_arr1 = x_initial(numDOF1:(numDOFScaffold-1));
    temp_arr2 = x_initial((almostAllVars+1):totVars);
    radius_squared = x_initial(1, 1) * x_initial(1, 1) + x_initial(numDOFScaffold, 1) * x_initial(numDOFScaffold, 1);
    if sum(temp_arr1 .* temp_arr1+temp_arr2 .* temp_arr2 > 0.8*radius_squared)
        break;
    end

    if norm(G)<gradTol
        break;
    end
    if (t*norm(d)<stepTol)
        break;
    end
    if (abs(E - Eprev) < Etol*(E + 1))
        if (EitersCounter >= numConvIters)
            break;
        else
            EitersCounter = EitersCounter + 1;
        end
    else
        EitersCounter = 0;
    end
end
arr1=[E];
GIF.newtonIters=iter;
GIF.NewtonFinalEnergy=E;
rate = abs(GIF.initTotalEnergy - E) / GIF.initTotalEnergy;

meanU = mean(x_initial(numDOF1:(numDOFScaffold-1)));
meanV = mean(x_initial((almostAllVars+1):totVars));
x_initial(numDOF1:(numDOFScaffold-1)) = x_initial(numDOF1:(numDOFScaffold-1)) - meanU;
x_initial((almostAllVars+1):totVars) = x_initial((almostAllVars+1):totVars) - meanV;
temp_arr1 = x_initial(numDOF1:(numDOFScaffold-1));
temp_arr2 = x_initial((almostAllVars+1):totVars);
maxRadiusSquared = max(temp_arr1 .* temp_arr1+temp_arr2 .* temp_arr2);
curOuterRadiusSquared = x_initial(1, 1) * x_initial(1, 1) + x_initial(numDOFScaffold, 1) * x_initial(numDOFScaffold, 1);
temp_scaling = 10 * sqrt(maxRadiusSquared / curOuterRadiusSquared);
% get full solution (with fixed)
solution=zeros(length(x_initial)+length(GIF.FixedIndices),1);
solution(setdiff(1:end,GIF.FixedIndices))=x_initial; % set free values
solution(GIF.FixedIndices)=GIF.FixedValues; % set fixed values
solution=[solution(1:n),solution(n+1:end)];
solution([1:numDOF1], :)=temp_scaling * solution([1:numDOF1], :);
a=solution([(numDOF1+1):numDOFScaffold], :);

GIF.UVonCones=a;

GIF.UVonCones1=[solution([1:numDOF1], :); circshift(solution([(numDOF1+1):numDOFScaffold], :),1,1)];
GIF.UVonCones2=[solution([1:numDOF1], :); flip(solution([(numDOF1+1):numDOFScaffold], :)); a];

isFound=0;
hole=0;
holeCoor=[0.0, 0.0];
for i=2:(length(a)-1)
    p1=complex(a(i-1, 1), a(i-1, 2));
    p2=complex(a(i, 1), a(i, 2));
    p3=complex(a(i+1, 1), a(i+1, 2));
    v1=p3-p2;
    v2=p1-p2;
    x1=real(v1);
    y1=imag(v1);
    x2=real(v2);
    y2=imag(v2);
    angle=atan2d(x1*y2-y1*x2,x1*x2+y1*y2);
    if(angle<179.0 && angle>0.0)
        hole=(p1+p2+p3)/3.0;
        isFound=1;
        break;
    end
end

if(isFound==1)
    holeCoor=[real(hole), imag(hole)];
else
    disp('Problem - cannot find a point inside the polygon');
end

GIF = rmfield(GIF,{'SF', 'SV', 'SJ_fz', 'SJ_fbz', 'SArea','TArea', 'EditedJ_fz', 'EditedJ_fbz','TJ_fz', 'TJ_fbz'});
clearvars -except GIF x_initial numR areaOriginal SPDHessian numDOF1 totVars almostAllVars iter holeCoor numFixedDOFScaffold timesRemeshed arr1 p gradTol1 stepTol1 Etol1 EitersCounter numConvIters terminationCondition rate lambda E LAST_EVAL_TO_COUT
