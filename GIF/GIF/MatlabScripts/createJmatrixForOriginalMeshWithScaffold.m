% create T_h - harmonic basis matrix for each triangle : dim: 3 X n X |F|
HarmonicBasisR1=full(HarmonicBasis(GIF.F(:,1),:)); % first row of T_h in each triangle dim: |F| X n 
HarmonicBasisR2=full(HarmonicBasis(GIF.F(:,2),:)); % second row of T_h in each triangle dim: |F| X n
HarmonicBasisR3=full(HarmonicBasis(GIF.F(:,3),:)); % third row of T_h in each triangle dim: |F| X n

% create T_fz and T_fbz - get the f_z of each triangle : dim: 1 X 3 X |F| (we use |F| X 3)

[tc,dblAc] = compute_perps(GIF.V, GIF.F', 3); %tc = |F| x 3 array in complex form, dblA =  array of doubled triangle areas
cf = -(1.0/2.0)./ dblAc;
t1r = tc(:,1); t2r =  tc(:,2); t3r =  tc(:,3);

J_fz=[cf.*conj(t1r) cf.*conj(t2r) cf.*conj(t3r)];% dim:|F| x 3
J_fbz=[cf.*t1r cf.*t2r cf.*t3r];% dim:|F| x 3

% create T_h_fz and T_h_fbz dim: 1 X n X |F| (we use |F| X n)
GIF.J_fz=bsxfun(@times, J_fz(:,1), HarmonicBasisR1)+bsxfun(@times, J_fz(:,2), HarmonicBasisR2)+bsxfun(@times, J_fz(:,3), HarmonicBasisR3);
GIF.J_fbz=bsxfun(@times, J_fbz(:,1), HarmonicBasisR1)+bsxfun(@times, J_fbz(:,2), HarmonicBasisR2)+bsxfun(@times, J_fbz(:,3), HarmonicBasisR3);

numR=size(GIF.J_fz, 1);
%get array of normalized triangles area
GIF.Area=dblAc/2;
internal_faces_indices = find(GIF.internalFacesIndicator == 1);
boundary_faces_indices = find(GIF.internalFacesIndicator == 0);
internal_area = GIF.Area(internal_faces_indices);
boundary_area = GIF.Area(boundary_faces_indices);
internal_area = internal_area ./ sum(internal_area);
boundary_area = boundary_area ./ sum(boundary_area);
boundary_area = 0.05 * boundary_area;
GIF.Area(internal_faces_indices) = internal_area;
GIF.Area(boundary_faces_indices) = boundary_area;
rate = 0;
lambda = 0;
areaOriginal=sum(GIF.Area);
iter=0;
numInternalTriangles = length(GIF.Area);
terminationCondition = true;