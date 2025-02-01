uvs = array1([1:GIF.sysVars], :);
uvs1 = uvs(:, 1)+1i*uvs(:, 2);
[k, triangleAreas, sigma_1, sigma_2] = extract_k(GIF.F1, GIF.V, uvs1);

sigma_1_squared = sigma_1 .^ 2;
sigma_2_squared = sigma_2 .^ 2;
part_a = (triangleAreas.' * (sigma_1_squared + sigma_2_squared));
part_b = (triangleAreas.' * (1 ./ sigma_1_squared + 1 ./ sigma_2_squared));
tot_area=sum(triangleAreas);
part_a1 = part_a/tot_area;
part_b1 = part_b/tot_area;
cur_energy = (1/2) * (part_a1+part_b1);

k_average = triangleAreas.' * k/tot_area;

disp(['symmetric Dirichlet energy is ',  num2str(cur_energy)]);
disp(['E_k is ', num2str(k_average)]);

energy_of_each_triangle = 0.5* (sigma_1_squared + sigma_2_squared + 1 ./ sigma_1_squared + 1 ./ sigma_2_squared);
energiesAndAreas = [energy_of_each_triangle, triangleAreas];
sortedEnergiesAndAreas = sortrows(energiesAndAreas, 1);
percent95length = round(0.95 * length(energiesAndAreas));
percent95triangles = sortedEnergiesAndAreas([1:percent95length], :);
totArea = sum(percent95triangles(:, 2));
percent95Energy = (percent95triangles(:, 1).' * percent95triangles(:, 2))/ totArea;

percent5triangles = sortedEnergiesAndAreas([percent95length+1:end], :);
totArea = sum(percent5triangles(:, 2));
percent5Energy = (percent5triangles(:, 1).' * percent5triangles(:, 2))/ totArea;

percent99length = round(0.99 * length(energiesAndAreas));
percent99triangles = sortedEnergiesAndAreas([1:percent99length], :);
totArea = sum(percent99triangles(:, 2));
percent99Energy = (percent99triangles(:, 1).' * percent99triangles(:, 2))/ totArea;

percent1triangles = sortedEnergiesAndAreas([percent99length+1:end], :);
totArea = sum(percent1triangles(:, 2));
percent1Energy = (percent1triangles(:, 1).' * percent1triangles(:, 2))/ totArea;

maxEnergy = max(energiesAndAreas(:, 1));
disp(['E_SD average over 95% of the triangles with the lowest distortion is ', num2str(percent95Energy)]);
disp(['E_SD average over 99% of the triangles with the lowest distortion is ', num2str(percent99Energy)]);
disp(['E_SD average over 5% of the triangles with the highest distortion is ', num2str(percent5Energy)]);
disp(['E_SD average over 1% of the triangles with the highest distortion is ', num2str(percent1Energy)]);
disp(['E_SD max is ', num2str(maxEnergy)]);

disp(['#Scaffold triangulations: ', num2str(timesRemeshed)]);
disp(['#Internal iterations: ', num2str(iter)]);

GIF.parametersMatrix(10,1) = percent99Energy;
GIF.parametersMatrix(9,1) = cur_energy;
GIF.parametersMatrix(11,1) = k_average;
GIF.parametersMatrix(16,1) = timesRemeshed;
GIF.parametersMatrix(17,1) = iter;
GIF.uvs = uvs;