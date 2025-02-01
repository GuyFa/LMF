function [  ] = saveObj( path, name, mesh )

obj.vertices = mesh.V;
obj.vertices_texture = mesh.uvs;
objects.data.vertices = mesh.F1;
objects.data.texture = mesh.F1;
obj.objects = objects;
obj.objects.type='f';

write_wobj(obj,[path '\\' name])

end

