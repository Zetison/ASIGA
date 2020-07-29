function plotModelInParaview(nurbs, options, plotMesh, appendName)
depricated('Use createParaviewFiles() instead')
extraXiPts = options.extraXiPts;
extraEtaPts = options.extraEtaPts;
extraZetaPts = options.extraZetaPts;

options = struct('name',[options.name, appendName]);

switch nurbs.type
    case '3Dsurface'
        options.celltype = 'VTK_QUAD';
        data = buildVisualizationMesh(nurbs, [extraXiPts, extraEtaPts]);
    case '3Dvolume'
        options.celltype = 'VTK_HEXAHEDRON';
        data = buildVisualizationMesh(nurbs, [extraXiPts, extraEtaPts, extraZetaPts]);
end

makeVTKfile(data, options);

if plotMesh
    Xi = nurbs.knots{1};
    Eta = nurbs.knots{2};
    noXiKnots = data.noXiKnots;
    noEtaKnots = data.noEtaKnots;
    XiVec = data.XiVec;
    EtaVec = data.EtaVec;
    createVTKmeshFiles(varCol, U, extraXiPts, extraEtaPts, extraZetaPts, options)
end
