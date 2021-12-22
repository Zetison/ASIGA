

% Define Lam? parameters
% lambda = nu*E/((1+nu)*(1-2*nu));
% mu = E/(2*(1+nu));
% 
% c_s = sqrt(E/((1-nu^2)*rho_s)); % Speed of sound in elastic material ??????????????????????????????????????
% c_1 = sqrt((lambda+2*mu)/rho_s);
% c_2 = sqrt(mu/rho_s);

alignWithAxis = 'Xaxis';

x_0 = [0; 0]; % The origin of the model
x_0wave = x_0; % Point at which the wave is directed (has purely real part here)
if useEllipse
    fact = 1.5;
else
    fact = 1;
end
fact2 = 2;
c_x = fact*fact2*R_o;
c_y = fact2*R_o;

f = sqrt(c_x^2-c_y^2);
[R_a, ~] = evaluateEllipticCoords(c_x,0,f);
% R_a = 10; %R_o + lambda/2;

rm = zeros(N+2,1);
for m = 1:N+2
    rm(m) = R_a + (m-1)*R_a;
end
% 
% stretchingFact = log(64/3)/log(2);
% stretchingFact = 10;

if i_mesh == 1
    noNewXiKnots = 0;
else
    noNewXiKnots = 1;
end
noNewXiKnots = ceil(pi/2*R_a)*(2^i_mesh-1);
noNewEtaKnots = ceil(R_a-R_o)*(2^i_mesh-1);
degreeElev = 0;
newMesh = 1;

if ~useEllipse

    fluid_o = getSphericalShellData_2D(R_o, R_a, alignWithAxis);
else
%     solid = getEllipseData(c_x/2,c_y/2);
%     fluid_o = embedSolidInEllipse(solid,c_x,c_y);
    solid = getEllipseData(1.5*R_o, R_o);
    fluid_o = embedSolidInEllipse(solid,c_x,c_y);
    
    
end

fluid_o = elevateNURBSdegree(fluid_o,[0 1]);
fluid_o = elevateNURBSdegree(fluid_o,[1 1]*degreeElev);
fluid_o = insertKnotsInNURBS(fluid_o,{insertUniform2(fluid_o.knots{1}, noNewXiKnots) ...
                                  insertUniform2(fluid_o.knots{2}, noNewEtaKnots)});


if useSolidDomain
    solid = getSphericalShellData(R, R_o, alignWithAxis);
    solid = elevateNURBSdegree(solid,[0 0 1]);
    solid = elevateNURBSdegree(solid,[1 1 1]*degreeElev);
    if i_mesh == 1
        noNewZetaKnots = 0;
    else
        noNewZetaKnots = 2;
    end
        
    solid = insertKnotsInNURBS(solid,{insertUniform2(solid.knots{1}, noNewXiKnots) ...
                                      insertUniform2(solid.knots{2}, noNewEtaKnots) ...
                                      insertUniform2(solid.knots{3}, noNewZetaKnots)});
end

if useInnerFluidDomain
    fluid_i = getSolidSphereData(R, alignWithAxis);
    fluid_i = elevateNURBSdegree(fluid_i,[0 0 1]);
    fluid_i = elevateNURBSdegree(fluid_i,[1 1 1]*degreeElev);
    if i_mesh == 1
        noNewZetaKnots = 0;
    else
        noNewZetaKnots = 2;
    end
        
    fluid_i = insertKnotsInNURBS(fluid_i,{insertUniform2(fluid_i.knots{1}, noNewXiKnots) ...
                                      insertUniform2(fluid_i.knots{2}, noNewEtaKnots) ...
                                      insertUniform2(fluid_i.knots{3}, noNewZetaKnots)});
%     error('this case is not implemented')
end


principalLength = 2*pi*R_o;

principalLengthXiDir = 2*pi*R_o;
principalLengthEtaDir = pi*R_o;




