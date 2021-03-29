function task = kirchApprTri(task,X)

noElems = task.varCol{1}.noElems;
element = task.varCol{1}.element;
tri = NaN(size(element,1),2,3);
P = task.varCol{1}.controlPts;
Eps = 1e2*eps;
for e = 1:noElems
    sctr = element(e,:);
    P1 = P(sctr(1),:);
    P2 = P(sctr(2),:);
    P3 = P(sctr(3),:);
    P4 = P(sctr(4),:);
    tri_e = NaN(1,2,3);
    if norm(P1-P2) < Eps
        tri_e(1,1,:) = element(e,[1,4,3]);
    elseif norm(P1-P3) < Eps
        tri_e(1,1,:) = element(e,[1,2,4]);
    elseif norm(P2-P4) < Eps || norm(P3-P4) < Eps
        tri_e(1,1,:) = element(e,[1,2,3]);
    else
        if norm(P2-P3) > norm(P1-P4)
            tri_e(1,1,:) = element(e,[1,2,4]);
            tri_e(1,2,:) = element(e,[1,4,3]);
        else
            tri_e(1,1,:) = element(e,[1,2,3]);
            tri_e(1,2,:) = element(e,[2,4,3]);
        end                                
    end
    tri(e,:,:) = tri_e;
end
tri = reshape(tri,size(tri,1)*size(tri,2),3);
tri(any(isnan(tri),2),:) = [];

%% Find h_max and store results
k = task.misc.omega/task.varCol{1}.c_f;
lambda = 2*pi./k;
task.varCol{1}.h_max = max([norm2(P(tri(:,1),:)-P(tri(:,2),:)); 
                    norm2(P(tri(:,1),:)-P(tri(:,3),:)); 
                    norm2(P(tri(:,2),:)-P(tri(:,3),:))]);
task.dofs = size(unique(tri,'rows','stable'),1);
task.varCol{1}.nepw = lambda./task.varCol{1}.h_max;
task.varCol{1}.noElems = size(tri,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% trisurf(tri,P(:,1),P(:,2),P(:,3), 'FaceColor', getColor(1))
% view(106,26) % sphere and cube
% axis off
% axis equal
% camlight
% ax = gca;               % get the current axis
% ax.Clipping = 'off';    % turn clipping off
% figureFullScreen(gcf)
% %                     
% export_fig(['../../graphics/sphericalShell/trianglesParm2_' num2str(task.varCol{1}.noElems)], '-png', '-transparent', '-r300')
%                     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



P_inc = task.misc.P_inc;
P1 = P(tri(:,1),:);
P2 = P(tri(:,2),:);
P3 = P(tri(:,3),:);

D0 = -2*1i*(P3*X.')*k;
D1 = -2*1i*((P1-P3)*X.')*k;
D2 = -2*1i*((P2-P3)*X.')*k;

g = (D1.*(1-exp(D2))+D2.*(exp(D1)-1))./D1./D2./(D1-D2);

indices1 = abs(D1) < Eps;
temp = D2(indices1);
g(indices1) = -(1+temp-exp(temp))./temp.^2;
indices2 = abs(D2) < Eps;
temp = D1(indices2);
g(indices2) = -(1+temp-exp(temp))./temp.^2;
indices3 = abs(D1-D2) < Eps;
temp = D2(indices3);
g(indices3) = (1-exp(temp)+temp.*exp(temp))./temp.^2;
g(and(indices1,indices2)) = 0.5;
normals = cross(P2-P1,P3-P1,2);
areas = norm2(normals); % (multiplied by 2)
n = normals./areas(:,[1,1,1]);
nX = n*X.';
areas = repmat(areas,1,size(nX,2));
areas(nX < 0) = 0;
task.ffp.p_h = -1i*k*P_inc/(2*pi).*sum(areas.*nX.*g.*exp(D0),1);

