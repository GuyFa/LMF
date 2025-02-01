function [k, triangleAreas, sigma_1, sigma_2] = extract_k(F, V, paramc)
%this function computes the conformal distortion for each triangle.
%larger than 1 value for k indicates that the triangle is flipped in the uv plane.
%it can accept as input a source mesh which is embedded in 3D but paramc is assumed to be complex (2D).
%the computation of the minimal and maximal singular values of the Jacobians is optional.
%the minimal singular value is signed. negative value indicates flipped triangle.

i_1 = F(:, 1);
i_2 = F(:, 2);
i_3 = F(:, 3);

p1 = V(i_1, :);
p2 = V(i_2, :);
p3 = V(i_3, :);

U1 = p3-p2;
U2 = p1-p2;

u1_length = sqrt(sum((U1.^2)'))';
dot = sum((U1.*U2)')';
cross_prod = cross(U1, U2, 2);

cross_prod = sqrt(sum((cross_prod.^2)'))';

triangleAreas = cross_prod/2;

e1 = u1_length;
e3 = -(dot + 1i*cross_prod)./u1_length;
e2 = -e1-e3;

t1 = -1i*e1;
t2 = -1i*e2;
t3 = -1i*e3;

f1 = paramc(i_1);
f2 = paramc(i_2);
f3 = paramc(i_3);

fzbarAbs = abs((f1.*t1 + f2.*t2 + f3.*t3) ./ (4*triangleAreas));
fzAbs = abs((f1.*conj(t1) + f2.*conj(t2) + f3.*conj(t3)) ./ (4*triangleAreas));

validIndices = (fzAbs ~= 0);

k(fzAbs == 0) = 1;
k(validIndices) = fzbarAbs(validIndices) ./ fzAbs(validIndices);
k = k';

if nargout > 2 
    sigma_1 = fzAbs + fzbarAbs;
    sigma_2 = fzAbs - fzbarAbs;
end

end