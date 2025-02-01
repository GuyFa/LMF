    areas = GIF.Area;
    areas = areas ./ sum(areas);
    
    [foldoversIndices, fz, fbz ,E]=symmetricDirichletEnergy(GIF.UVonCones,GIF.J_fz,GIF.J_fbz,areas,GIF.FixedIndices,GIF.FixedValues,true, 1);
    GIF.initTotalEnergy=E;
    GIF.initEnergy=E;
    % optimize scale
    scale = optimizeSymmetricDirichletEnergyByGlobalScaling( fz, fbz, GIF.Area );
    GIF.UVonCones=scale*GIF.UVonCones;
    GIF.FixedValues=scale*GIF.FixedValues;
    
    
    [~, ~, ~, E]=symmetricDirichletEnergy(GIF.UVonCones,GIF.J_fz,GIF.J_fbz,areas,GIF.FixedIndices,GIF.FixedValues,true, 1);
    x9=zeros(length(GIF.UVonCones)+length(GIF.FixedIndices),1);
    x9(setdiff(1:end,GIF.FixedIndices))=GIF.UVonCones;
    x9(GIF.FixedIndices)=GIF.FixedValues;
    x10=[x9(1:(length(x9)/2)), x9(((length(x9)/2)+1):length(x9))];
