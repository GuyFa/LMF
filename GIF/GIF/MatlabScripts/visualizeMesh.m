function [] = visualizeMesh(mesh)
figure;
trisurf(mesh.F1,mesh.V(:,1),mesh.V(:,2),mesh.V(:,3));
hold on
axis equal
grid off
axis off
axis vis3d

plot3(mesh.borderEdges(:,1), mesh.borderEdges(:,2), mesh.borderEdges(:,3),'Color','r','LineWidth',3);
end

