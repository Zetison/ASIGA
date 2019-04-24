
%% Set optimal parameters for radial distributions of knots (from rigid sphere analysis)
c1 = 0.009435630480941;
c2 = 0.658633686954660;
zeta_u4 = 0.416774132047842;
R_a = 6;
% R_a = 1.2*R_o;

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


Upsilon = 0;
% R_a = R_o + lambda/2;

rm = zeros(N+2,1);
for m = 1:N+2
    rm(m) = R_a + (m-1)*R_a;
end
% 
% stretchingFact = log(64/3)/log(2);
% stretchingFact = 10;

noNewXiKnots = 2^(i_mesh-1)-1; % 8*i_mesh
noNewEtaKnots = 2^(i_mesh-1)-1;
noNewZetaKnots = 2^(i_mesh-1)-1;

% t = R_o*(1.2-1);
% f_arc = @(s) sqrt(R_o^2*sin(s).^2+R_o^2*cos(s).^2);
% noNewXiKnots = 3*(2^(i_mesh-1)-1);
% noNewEtaKnots = 3*round(integral(f_arc,0,pi/2)/(R_o*pi/2)*2^(i_mesh-1));
% noNewZetaKnots = 2*round(t/(R_o*pi/2)*(2^(i_mesh-1)+1));

fluid_o = getSphericalShellData(R_o, R_a, alignWithAxis);

fluid_o = elevateNURBSdegree(fluid_o,[0 0 1]);
fluid_o = elevateNURBSdegree(fluid_o,[1 1 1]*degreeElev);

% radZknots = insertUniform2(fluid_o.knots{3}, noNewZetaKnots);

zeta_vec = linspace2(0,1,noNewZetaKnots);
if true
    coeffs = [0, c1, c2, 1];
    Zeta_u = [0 0 0 zeta_u4 1 1 1];
else
    coeffs = [0, 0.5, 1];    
    Zeta_u = [0 0 0.5 1 1];
end

newZetaKnots = ZetaKnotsDistribution(coeffs, Zeta_u, zeta_vec);

fluid_o = insertKnotsInNURBS(fluid_o,{insertUniform2(fluid_o.knots{1}, noNewXiKnots) ...
                                  insertUniform2(fluid_o.knots{2}, noNewEtaKnots) ...
                                  newZetaKnots});

if useSolidDomain
    solid = getSphericalShellData(R_i, R_o, alignWithAxis);
    solid = elevateNURBSdegree(solid,[0 0 1]);
    solid = elevateNURBSdegree(solid,[1 1 1]*degreeElev);
        
    solid = insertKnotsInNURBS(solid,{insertUniform2(solid.knots{1}, noNewXiKnots) ...
                                      insertUniform2(solid.knots{2}, noNewEtaKnots) ...
                                      insertUniform2(solid.knots{3}, noNewZetaKnots)});
end

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


