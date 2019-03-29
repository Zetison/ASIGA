function min_d = minDistNURBS(x,nurbsCol)
tic
if ~iscell(nurbsCol)
    nurbsCol = {nurbsCol};
end
min_d = inf;
for i_nurbs = 1:numel(nurbsCol)
    nurbs = nurbsCol{i_nurbs};
    dXdxi = @(pt) evaluateNURBS_deriv2(nurbs, pt.', 'xi');
    dXdeta = @(pt) evaluateNURBS_deriv2(nurbs, pt.', 'eta');
    d2Xdxi2 = @(pt) evaluateNURBS_2ndDeriv2(nurbs, pt.', 'xi');
    d2Xdeta2 = @(pt) evaluateNURBS_2ndDeriv2(nurbs, pt.', 'eta');
    d2Xdxideta = @(pt) evaluateNURBS_2ndDeriv2(nurbs, pt.', 'xieta');
    X = @(pt) evaluateNURBS(nurbs, pt);
    f = @(pt) [dot((x-X(pt)),dXdxi(pt));
               dot((x-X(pt)),dXdeta(pt))];

    repKnotsXi = repKnots(nurbs.degree(1),nurbs.knots{1}(2:end-1));
    repKnotsEta = repKnots(nurbs.degree(2),nurbs.knots{2}(2:end-1));

    bnd = zeros(2,2,(numel(repKnotsXi)-1)*(numel(repKnotsEta)-1));
    counter = 1;
    for j = 1:numel(repKnotsEta)-1
        for i = 1:numel(repKnotsXi)-1
            bnd(1,:,counter) = [repKnotsXi(i),repKnotsXi(i+1)];
            bnd(2,:,counter) = [repKnotsEta(j),repKnotsEta(j+1)];
            counter = counter + 1;
        end
    end

    nptsXi = 2;
    nptsEta = 1;
    extremas = zeros(size(bnd,3)*nptsXi*nptsEta,2);
    noItrs = 100;
    counter = 1;
    for i = 1:size(bnd,3)
        xi_arr = linspace2(bnd(1,1,i),bnd(1,2,i),nptsXi);
        eta_arr = linspace2(bnd(2,1,i),bnd(2,2,i),nptsEta);
        for j = 1:nptsEta
            for l = 1:nptsXi
                [extremas(counter,:),nrItr] = newtonsMethodND(f,@(pt)jac(pt,x),[xi_arr(l),eta_arr(j)].',noItrs,1e-14,bnd(:,:,i));
    %             nrItr
                counter = counter + 1;
            end
        end
    end

    d = zeros(size(extremas,1),1);
    for i = 1:size(extremas,1)
        d(i) = norm(x-X(extremas(i,:)));
    end
    [d,I] = min(d);
    if d < min_d
        min_d = d;
    end
%     pt = extremas(I,:);

end
% disp(['One evaluation took ', num2str(toc), ' seconds!']);

function J = jac(pt,x)

J = zeros(2);
J(1,1) = -norm(dXdxi(pt))^2 + dot((x-X(pt)),d2Xdxi2(pt));
J(1,2) = -dot(dXdxi(pt),dXdeta(pt)) + dot((x-X(pt)),d2Xdxideta(pt));
J(1,2) = J(2,1);
J(2,2) = -norm(dXdeta(pt))^2 + dot((x-X(pt)),d2Xdeta2(pt));
end
end

function repKnots = repKnots(p,Xi)
if p == 1
    repKnots = Xi;
    return
end
i = 1;
repKnots = zeros(size(Xi));
counter = 0;
while i < numel(Xi)
    xi = Xi(i);
    m = sum(Xi == xi);
    if m == p
        counter = counter + 1;
        repKnots(counter) = xi;
    end
    i = i + m;
end
repKnots = repKnots(1:counter);
end