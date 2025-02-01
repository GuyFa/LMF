
% create T_h - harmonic basis matrix for each triangle : dim: 3 X n X |F|
ScaffoldHarmonicBasisR1=full(ScaffoldHarmonicBasis(GIF.SF(:,1),:)); % first row of T_h in each triangle dim: |F| X n 
ScaffoldHarmonicBasisR2=full(ScaffoldHarmonicBasis(GIF.SF(:,2),:)); % second row of T_h in each triangle dim: |F| X n
ScaffoldHarmonicBasisR3=full(ScaffoldHarmonicBasis(GIF.SF(:,3),:)); % third row of T_h in each triangle dim: |F| X n
% create T_fz and T_fbz - get the f_z of each triangle : dim: 1 X 3 X |F| (we use |F| X 3)

[tc,dblAc] = compute_perps(GIF.SV, GIF.SF', 3); %tc = |F| x 3 array in complex form, dblA =  array of doubled triangle areas
cf = -(1.0/2.0)./ dblAc;
t1r = tc(:,1); t2r =  tc(:,2); t3r =  tc(:,3);

J_fz=[cf.*conj(t1r) cf.*conj(t2r) cf.*conj(t3r)];% dim:|F| x 3
J_fbz=[cf.*t1r cf.*t2r cf.*t3r];% dim:|F| x 3

% create T_h_fz and T_h_fbz dim: 1 X n X |F| (we use |F| X n)
GIF.SJ_fz=bsxfun(@times, J_fz(:,1), ScaffoldHarmonicBasisR1)+bsxfun(@times, J_fz(:,2), ScaffoldHarmonicBasisR2)+bsxfun(@times, J_fz(:,3), ScaffoldHarmonicBasisR3);
GIF.SJ_fbz=bsxfun(@times, J_fbz(:,1), ScaffoldHarmonicBasisR1)+bsxfun(@times, J_fbz(:,2), ScaffoldHarmonicBasisR2)+bsxfun(@times, J_fbz(:,3), ScaffoldHarmonicBasisR3);


numDOF =size(GIF.J_fz, 2);%82
numDOFScaffold=size(GIF.SJ_fz, 2);%410

numFixedDOFScaffold=numDOFScaffold-numDOF;%328
GIF.EditedJ_fz=[zeros(numR, numFixedDOFScaffold), GIF.J_fz];
GIF.EditedJ_fbz=[zeros(numR, numFixedDOFScaffold), GIF.J_fbz];
GIF.SJ_fz=[GIF.SJ_fz(:, [1:numFixedDOFScaffold]), flip(GIF.SJ_fz(:, [(numFixedDOFScaffold+1):numDOFScaffold]), 2)];
GIF.SJ_fbz=[GIF.SJ_fbz(:, [1:numFixedDOFScaffold]), flip(GIF.SJ_fbz(:, [(numFixedDOFScaffold+1):numDOFScaffold]), 2)];

if iter == 0
    lambda = E/200;
else
    lambda = max(1e-9, rate);
end
GIF.SArea = dblAc/2;
GIF.SArea = GIF.SArea ./ sum(GIF.SArea);
GIF.SArea = lambda * GIF.SArea;
GIF.TArea=[GIF.Area;GIF.SArea];
GIF.TArea = GIF.TArea ./ sum(GIF.TArea);
GIF.TJ_fz=[GIF.EditedJ_fz;GIF.SJ_fz];
GIF.TJ_fbz=[GIF.EditedJ_fbz;GIF.SJ_fbz];
numDOF1=numDOFScaffold-numDOF;

%  clear ScaffoldHarmonicBasis ScaffoldHarmonicBasisR1 ScaffoldHarmonicBasisR2 ScaffoldHarmonicBasisR3 tc dblAc cf t1r t2r t3r J_fz J_fbz