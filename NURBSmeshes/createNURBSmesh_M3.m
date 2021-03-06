function varCol = createNURBSmesh_M3(varCol, M, degree)
if isfield(varCol{1}, 'R1')
    R1 = varCol{1}.R1;
    R2 = varCol{1}.R2;
else
    R1 = varCol{1}.R;
    R2 = varCol{1}.R;
end
Xi = [0,0,0,1,1,2,2,3,3,3]/3;
t = varCol{1}.t;
L = varCol{1}.L;
parm = varCol{1}.parm;
x_0 = [-L/2-(R2-R1)/2, 0, 0]; % The origin of the model
alignWithAxis = 'Xaxis';
switch varCol{1}.method
    case {'IE','IENSG','ABC','MFS'}
        % A_2 maps the x-axis to the z-axis, the y-axis to the x-axis, and
        % the z-axis to the y-axis
        varCol{1}.x_0 = x_0;
        A_2 = [0 1 0;
               0 0 1;
               1 0 0];
        varCol{1}.A_2 = A_2;
        varCol{1}.alignWithAxis = alignWithAxis;
end

if varCol{1}.boundaryMethod
    c_z = (L+R1+R2)/2;
    c_xy = sqrt(R1*R2*(L^2 - (R2 - R1)^2))/L;
    varCol{1}.c_z = c_z;
    varCol{1}.c_xy = c_xy;

    Upsilon = sqrt(c_z^2-c_xy^2);

    chimin = 24.4;
    chimax = 25.7;

    solid = getBeTSSiM3Data('R1', R1, 'R2', R2, 't', t, 'L', L, 'parm', parm, 'Xi', Xi);
    if varCol{1}.refineThetaOnly
        if parm ~= 1
            error('Must have parm = 1 for pure theta refinement')
        end
        solid = makeUniformNURBSDegree(solid,[2,degree,degree]);
        solid = insertKnotsInNURBS(solid,[0,2^(M-1)-1,0]);
    else
        solid = makeUniformNURBSDegree(solid,degree);
        Imap{1} = [R2*pi/2,R1*pi/2,(R2-t)*pi/2,(R1-t)*pi/2];
        solid = refineNURBSevenly(solid,(2^(M-1)-1)/(R2*pi/2),Imap);
    end

    fluid = subNURBS(solid, 'at', [0,0;0,0;0,1]);
    varCol{1}.patchTop = getPatchTopology(fluid);
    if numel(varCol) > 2
        fluid_i = subNURBS(solid, 'at', [0,0;0,0;1,0]);
        varCol{3}.patchTop = getPatchTopology(fluid_i);
    end
    varCol_dummy.dimension = 1;
    varCol_dummy.nurbs = fluid;
    varCol_dummy = findDofsToRemove(generateIGAmesh(convertNURBS(varCol_dummy)));
    
    varCol{1}.elemsOuter = 1:varCol_dummy.noElems;
    varCol{1}.noDofsInner = 0;
    varCol{1}.noElemsInner = 0;
else
    if strcmp(varCol{1}.model,'MS')
        Xi = [0,0,0,1,1,2,2,3,3,4,4,4]/4;
        c_z = 1.2*(L+2*R1)/2; % 30
        c_xy = 1.3*R1; % 2.5, 3.75
        Upsilon = sqrt(c_z^2-c_xy^2);
        eta1 = 0.37;
        eta2 = 1-eta1;
        noNewZetaKnots = max(2^(M-1)/4-1,0);
        nn = 2^(M-1)-1;
    else
        R_max = max(R1,R2);
        s = 0.25 + 0.05*(L-R_max)/R_max;
        c_z = (L+R1+R2)/2 + s*R_max;
        c_xy = R_max + s*R_max;
        Upsilon = sqrt(c_z^2-c_xy^2);
        [~, theta1, ~] = evaluateProlateCoords(([0,0,R1]-x_0)*A_2.',Upsilon);
        [~, theta2, ~] = evaluateProlateCoords(([-L,0,R1]-x_0)*A_2.',Upsilon);
        eta1 = theta1/pi;
        eta2 = theta2/pi;
    end

    varCol{1}.r_a = evaluateProlateCoords([0,0,c_z],Upsilon);
    
    chimin = NaN;
    chimax = NaN;

    solid = getBeTSSiM3Data('R1', R1, 'R2', R2, 't', t, 'L', L, 'parm', parm, 'Xi', Xi);
    solid = [glueNURBS(solid(1:3),1),glueNURBS(solid(4:6),1),glueNURBS(solid(7:9),1)];
    if varCol{1}.refineThetaOnly
        if parm ~= 1
            error('Must have parm = 1 for pure theta refinement')
        end
        degreeVec = [2,degree,degree];
        refIndices = 2:3;
    else
        degreeVec = degree;
        refIndices = 1:3;
    end
    
    ellipsoid = getEllipsoidData('C',[c_z,c_xy,c_xy],'alignWithAxis',alignWithAxis,'x_0',x_0, 'alpha', 0, ...
                                 'Xi', Xi, 'Eta', [0,0,0,eta1,eta1,eta2,eta2,1,1,1]);
	
    if strcmp(varCol{1}.model,'MS')
        fluid = loftNURBS({subNURBS(solid,'at',[0,0;0,0;0,1]),explodeNURBS(ellipsoid,2)});
        fluid = makeUniformNURBSDegree(fluid,degree);
        fluid = glueNURBS([glueNURBS(fluid(1:4),1),glueNURBS(fluid(5:8),1),glueNURBS(fluid(9:12),1)],2);
        fluid = insertKnotsInNURBS(fluid,[nn,nn,noNewZetaKnots]);
    else
        ellipsoid = makeUniformNURBSDegree(ellipsoid,degreeVec(1:2));

        solid = makeUniformNURBSDegree(solid,degreeVec);
%         nin = [solid{1}.number(2),solid{4}.number(2),solid{7}.number(2)]-(solid{1}.degree(2)+1);
%         nin = round([R1*pi/2,L,R2*pi/2]*(2^(M-1)-1)/(R1*pi/2));
%         etas = [0,eta1,eta2,1];
%         newKnots = cell(1,3);
%         for i = 1:numel(nin)
%             Ltot = NURBSarcLength(ellipsoid{1},etas(i),etas(i+1),[0,NaN],2);
%             deta = (etas(i+1)-etas(i))/(nin(i)+1);
%             newKnots{i} = zeros(1,nin(i));
%             for j = 1:nin(i)
%                 if j == 1
%                     prevEta = etas(i);
%                 else
%                     prevEta = newKnots{i}(j-1);
%                 end
%                 L_j = @(eta) NURBSarcLength(ellipsoid{1},prevEta,eta,[0,NaN],2);
%                 f_j = @(eta) L_j(eta)-Ltot/(nin(i)+1);
%                 dfdeta_j = @(eta) dfdeta_(ellipsoid{1},[0,eta]);
%                 newKnots{i}(j) = newtonsMethod(f_j,dfdeta_j,prevEta+deta,100,eps,[prevEta,etas(i+1)]);
%             end
%             solid(i) = insertKnotsInNURBS(solid(i),{[], (newKnots{i}-etas(i))/(etas(i+1)-etas(i)), []});
%         end
%         ellipsoid = insertKnotsInNURBS(ellipsoid,{[],[newKnots{1},newKnots{2},newKnots{3}]});
%         
        fluid = loftNURBS({subNURBS(solid,'at',[0,0;0,0;0,1]),explodeNURBS(ellipsoid,2)});
        
        Imap{1} = [R2*pi/2, 10.896638128449496];
        Imap{2} = [4.33705, 0.813387090943184];
        Imap{3} = [L, 68.942341105179977];
        fluid = refineNURBSevenly(fluid,(2^(M-1)-1)/(R2*pi/2),Imap,0,2:3,0);
        solid = refineNURBSevenly(solid,(2^(M-1)-1)/(R2*pi/2),{},0,2:3,0); 
        uniqueEta = unique(fluid{1}.knots{2});
        extraEtaKnot = uniqueEta(2)/2;
        fluid(1) = insertKnotsInNURBS(fluid(1),{[],1-extraEtaKnot,[]});
        fluid(3) = insertKnotsInNURBS(fluid(3),{[],extraEtaKnot,[]});
        solid(1) = insertKnotsInNURBS(solid(1),{[],1-extraEtaKnot,[]});
        solid(3) = insertKnotsInNURBS(solid(3),{[],extraEtaKnot,[]});
%         ellipsoid = subNURBS(fluid,'at',[0,0;0,0;0,1]);
%         plotNURBS(ellipsoid,'colorFun',@(x) log10(abs(sum((x-x_0).^2./[c_z,c_xy,c_xy].^2)-1)))
    end
end
if parm == 2 && degree < 4
    warning(['parm=2 requires degree >= 4. Using degree=4 instead of degree=' num2str(degree)])
end
varCol{1}.chimin = chimin;
varCol{1}.chimax = chimax;
varCol{1}.L_gamma = L + R1 + R2;
varCol{1}.Upsilon = Upsilon;

varCol{1}.nurbs = fluid;
if numel(varCol) > 1
    varCol{2}.nurbs = solid;
end
if numel(varCol) > 2
    varCol{3}.nurbs = fluid_i;
end

function dfdeta_ = dfdeta_(nurbs,xi)
X = cell(1,3);
[X{:}] = evaluateNURBS(nurbs, xi, 1);
dfdeta_ = norm(X{3});