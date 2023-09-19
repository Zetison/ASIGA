function nurbs = convexifyNURBS(nurbs)

controlPts = [];
if nurbs{1}.d_p ~= 2
    error('Case not implemented')
end
options.Eps = 1e-10;
% options.avg_v_n_threshholdAngle = Inf; % Allways use the average normal vectors

nurbs = cleanNURBS(nurbs,[],options.Eps); % Important for convhull as it is sensitive to eps differences in points that should be considered the same

S2Vobj = surfaceToVolume(nurbs,options);
for patch = 1:numel(S2Vobj.nurbs_bdry)
    controlPts = [controlPts; S2Vobj.nurbs_bdry{patch}.coeffs(1:3,:).'];
end
TRI = convhull(controlPts,'simplify',true);
P1 = controlPts(TRI(:,1),:);
P2 = controlPts(TRI(:,2),:);
P3 = controlPts(TRI(:,3),:);
normalsTRI = normalize_lp(cross(P2-P1,P3-P1).').';
v1 = P2-P1;
v1s = sum(v1.^2,2);
v2 = P3-P1;
v2s = sum(v2.^2,2);
v1dv2 = dot(v1,v2,2);
D = v1s.*v2s-v1dv2.^2;

for patch = 1:numel(S2Vobj.nurbs_bdry)
    face = faceFromNormals(S2Vobj,patch);
    normals = face{1}.coeffs(1:3,:).';
    X = S2Vobj.nurbs_bdry{patch}.coeffs(1:3,:).';
    
    for i = 1:size(X,1)
        d = dot(X(i,:)-P1,normalsTRI,2);
        nnTRI = sum(normals(i,:).*normalsTRI,2);
        d(nnTRI < -options.Eps) = -Inf; % These points are on the oposite side of the closed surface
        p = X(i,:) - normalsTRI.*d;

        pmP1 = p-P1;
        pmP1dv1 = dot(pmP1,v1,2);
        pmP1dv2 = dot(pmP1,v2,2);
        xi2 = (v2s.*pmP1dv1-v1dv2.*pmP1dv2)./D;
        xi3 = (v1s.*pmP1dv2-v1dv2.*pmP1dv1)./D;
        hit = and(and(xi2 >= -options.Eps, xi3 >= -options.Eps), xi2+xi3 < 1+options.Eps);
        [~,min_d_I] = max(d(hit));
%         [~,min_d_I] = max(d(hit)./abs(nnTRI(hit))); % Scale in order to encourage aligned normal vectors
        p_temp = p(hit,:);
        if isempty(p_temp)
            warning('Could not find a proper intersection point')
        else
            S2Vobj.nurbs_bdry{patch}.coeffs(1:3,i) = p_temp(min_d_I,:);
        end
    end
end
nurbs = S2Vobj.nurbs_bdry;
