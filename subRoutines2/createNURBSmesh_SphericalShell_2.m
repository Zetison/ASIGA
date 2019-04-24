if useSolidDomain
    % Define Lam? parameters
    lambda = nu*E/((1+nu)*(1-2*nu));
    mu = E/(2*(1+nu));

    c_s = sqrt(E/((1-nu^2)*rho_s)); % Speed of sound in elastic material ??????????????????????????????????????
    c_1 = sqrt((lambda+2*mu)/rho_s);
    c_2 = sqrt(mu/rho_s);
end

alignWithAxis = 'Zaxis';

x_0 = [0; 0; 0]; % The origin of the model
x_0wave = x_0; % Point at which the wave is directed (has purely real part here)


Upsilon = 0; %0.7*R_o;
chimin = 0.9*R_o;
chimax = 2.1*R_o;
c_z = R_o;
c_xy = sqrt(c_z^2-Upsilon^2);

noNewXiKnots = 3*(2^(i_mesh-1)-1); % 8*i_mesh
noNewEtaKnots = 3*(2^(i_mesh-1)-1);
noNewZetaKnots = 2^(i_mesh-1)-1;

solid = getSphericalShellData(R_i, R_o, alignWithAxis);
solid = elevateNURBSdegree(solid,[0 0 1]);
solid = elevateNURBSdegree(solid,[1 1 1]*degreeElev);

solid = insertKnotsInNURBS(solid,{insertUniform2(solid.knots{1}, noNewXiKnots) ...
                                  insertUniform2(solid.knots{2}, noNewEtaKnots) ...
                                  insertUniform2(solid.knots{3}, noNewZetaKnots)});


if useInnerFluidDomain
    fluid_i = getSolidSphereData(R_i, alignWithAxis);
    fluid_i = elevateNURBSdegree(fluid_i,[0 0 1]);
    fluid_i = elevateNURBSdegree(fluid_i,[1 1 1]*degreeElev);
        
    fluid_i = insertKnotsInNURBS(fluid_i,{insertUniform2(fluid_i.knots{1}, noNewXiKnots) ...
                                      insertUniform2(fluid_i.knots{2}, noNewEtaKnots) ...
                                      insertUniform2(fluid_i.knots{3}, noNewZetaKnots)});
%     error('this case is not implemented')
end


principalLength = 2*pi*R_o;

principalLengthXiDir = 2*pi*R_o;
principalLengthEtaDir = pi*R_o;


L_gamma = 2*R_o;


