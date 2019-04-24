function [tri, X, pts, pts2] = triangulateCrossSection(nurbs, parm_pt, dir, delta, xb, yb, zb)


if numel(zb) == 1
    pplane = 3;
elseif numel(yb) == 1
    pplane = 2;
elseif numel(xb) == 1
    pplane = 1;
end

totalArcLength = arcLength(nurbs,0,1,parm_pt,dir);
nIntervals = floor(totalArcLength/delta);
 
pts = zeros(1,nIntervals+1);
pt_prev = 0;
i = 1;
while 1 - pts(i) > 1e-10
    fun = @(pt) abs(arcLength(nurbs,pt_prev,pt,parm_pt,dir) - delta);
    pts(i+1) = fminsearchbnd(fun,pt_prev,pt_prev,1);
    pt_prev = pts(i+1);
    i = i + 1;
end
pts = sort(unique([nurbs.knots{dir},pts]));
pts(end) = [];
v = zeros(length(pts),3);
for i = 1:length(pts)
    switch dir
        case 1
            v(i,:) = evaluateNURBS(nurbs, [pts(i), parm_pt(2)]);
        case 2
            v(i,:) = evaluateNURBS(nurbs, [parm_pt(1), pts(i)]);
    end
end
if dir == 2    
    pts2 = [1, fliplr(pts(2:end))];
    v2 = zeros(length(pts2),3);
    for i = 1:length(pts2)
        v2(i,:) = evaluateNURBS(nurbs, [parm_pt(2), pts2(i)]);
    end
    v = [v; v2];
else
    pts2 = [];
end

switch pplane
    case 3
        p1 = v(:,1:2);
        p2 = [p1(2:end,:); p1(1,:)];
        p3 = p2-p1;
        p4 = [-p3(:,2), p3(:,1)];
        
        v2 = (p1+p2)/2 + sqrt(3)/2*delta*p4./repmat(norm2(p4),1,2);
        X = triangulateRectangle([xb(1), yb(1)],[xb(2), yb(2)],round((xb(2)-xb(1))/delta));
        
        in = inpolygon(X(:,1),X(:,2),v2(:,1),v2(:,2));
        X(in,:) = [];
        X = [v(:,1),v(:,2); v2; X];
    case 2
        p1 = v(:,[1 3]);
        p2 = [p1(2:end,:); p1(1,:)];
        p3 = p2-p1;
        p4 = [p3(:,2), -p3(:,1)];
        
        v2 = (p1+p2)/2 + sqrt(3)/2*delta*p4./repmat(norm2(p4),1,2);
        X = triangulateRectangle([xb(1), zb(1)],[xb(2), zb(2)],round((xb(2)-xb(1))/delta));
        in = inpolygon(X(:,1),X(:,2),v2(:,1),v2(:,2));
        X(in,:) = [];
        X = [v(:,1),v(:,3); v2; X];
    case 1
        p1 = v(:,2:3);
        p2 = [p1(2:end,:); p1(1,:)];
        p3 = p2-p1;
        p4 = [p3(:,2), -p3(:,1)];
        v2 = (p1+p2)/2 + sqrt(3)/2*delta*p4./repmat(norm2(p4),1,2);
        
        X = triangulateRectangle([yb(1), zb(1)],[yb(2), zb(2)],round((yb(2)-yb(1))/delta));
        in = inpolygon(X(:,1),X(:,2),v2(:,1),v2(:,2));
        X(in,:) = [];
        X = [v(:,2),v(:,3); v2; X];
end
nkpts = size(v,1);
C = [1:nkpts;
     2:nkpts, 1].';
 
tri = delaunayTriangulation(X(:,1),X(:,2), C);
trianglesInside = isInterior(tri);
tri = tri(~trianglesInside,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% keyboard
% close all
% figure(1)
% % triplot(tri(~IO, :),tri.Points(:,1), tri.Points(:,2),'LineWidth', 2)
% triplot(tri,X(:,1),X(:,2))
% % for i = 1:size(X,1)
% %     text(X(i,1),X(i,2),num2str(i), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')
% % end
% % for i = 1:size(tri,1)
% %     x = mean(X(tri(i,:),:));
% %     text(x(1),x(2),num2str(i), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle','color','red')
% % end
% axis equal
% hold on
% 
% % plot(v(:,2),v(:,3),'*','color','green')
% % plot(X1(:,1),X1(:,2),'o','color','black')
% % plot(X(indices,1),X(indices,2),'*','color','green')
% % 
% switch pplane
%     case 3
%         plot(v(:,1),v(:,2),'*-','color','red')
%         plot(v2(:,1),v2(:,2),'*-','color','green')
%     case 2
%         plot(v(:,1),v(:,3),'*-','color','red')
%         plot(v2(:,1),v2(:,2),'*-','color','green')
%     case 1
%         plot(v(:,2),v(:,3),'*-','color','red')
%         plot(v2(:,1),v2(:,2),'*-','color','green')
% end
% 
% keyboard
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = evaluateNURBS(nurbs,parm_pt);
switch pplane
    case 3
        X = [X, ones(size(X,1),1)*x(3)];
    case 2
        X = [X(:,1), ones(size(X,1),1)*x(2), X(:,2)];
    case 1
        X = [ones(size(X,1),1)*x(1), X];
end

 
function I = arcLength(nurbs,a,b,parm_pt,dir)
 
switch dir
    case 1
        fun = @(pt) integrand(nurbs,dir,[pt, parm_pt(2)]);
    case 2
        fun = @(pt) integrand(nurbs,dir,[parm_pt(1), pt]);
end
 
n = 20;
[W1D,Q1D] = gaussianQuadNURBS(n); 
I = 0;
for i = 1:n
    I = I + fun((b-a)/2*(1+Q1D(i))+a)*W1D(i);
end
I = I*(b-a)/2;

 
function I = integrand(nurbs,dir,parm_pt)
 
 
switch dir
    case 1
        [~, I] = evaluateNURBS_deriv(nurbs, parm_pt);
    case 2
        [~, ~, I] = evaluateNURBS_deriv(nurbs, parm_pt);
end
 
I = norm(I);