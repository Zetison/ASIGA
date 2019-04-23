
if varCol.boundaryMethod
    alignWithAxis = 'Zaxis';
    c_z = 0.98*(L+R_o)/2;
    c_xy = 0.99*R_o/2; % 2.5, 3.75; 
    Upsilon = sqrt(c_z^2-c_xy^2);

    chimin = 21.07;
    chimax = 23.1;

    eta2 = 0.8*(c_z/30)^(1-1.0);
    eta1 = 0.265*(c_z/30)^(1-1.6);


    solid = getModel1Data(R_i,R_o,L, x_0, alignWithAxis, eta1, eta2);
    solid = elevateNURBSdegree(solid,[0 0 1]);
    solid = elevateNURBSdegree(solid,[1 1 1]*degreeElev);

    noNewXiKnots = 2^(M-1)-1; % 8*i_mesh
    noNewEtaKnotsEnds = noNewXiKnots;
    noNewEtaKnotsMiddle = floor(noNewXiKnots*L/((R_o + R_o*pi/2)/2));

    noNewZetaKnots = 2^(M-1)-1;
    solid = insertKnotsInNURBS(solid,{insertUniform2(solid.knots{1}, noNewXiKnots) ...
                                      [insertUniform2([0 eta1], noNewEtaKnotsEnds); ...
                                       insertUniform2([eta1 eta2], noNewEtaKnotsMiddle); ...
                                       insertUniform2([eta2 1], noNewEtaKnotsEnds)] ...
                                      insertUniform2(solid.knots{3}, noNewZetaKnots)});

    principalLengthXiDir = pi*R_o+2*R_o;
    principalLengthEtaDir = R_o*pi/2 + R_o + L;

    fluid = extractOuterSurface(solid);
    
    varCol.patchTop{1} = [ones(4,1),zeros(4,1)];
    varCol.patchTop{1}(2,2) = NaN;
    varCol.patchTop{1}(4,2) = NaN;
end
switch method
    case {'IE','IENSG'}
        varCol.x_0 = [-R_o/2, 0, (L-R_o)/2]; % The origin of the model
        % A_2 maps the x-axis to the z-axis, the y-axis to the x-axis, and
        % the z-axis to the y-axis
        varCol.A_2 = [0 1 0;
                      0 0 1;
                      1 0 0];
end