function [] = visualizeUVLayout(mesh)
figure;
triplot(mesh.F1, mesh.uvs(:, 1), mesh.uvs(:, 2), 'LineStyle', '-', 'Color', 'k','LineWidth', 0.1);
hold on
axis equal
grid off
axis off
axis vis3d

end


