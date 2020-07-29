function d = minDist(v,X,dXdxi,dXdeta,d2Xdxi2,d2Xdeta2,d2Xdxideta,bnd)
tic

f = @(pt) [dot((v-X(pt)),dXdxi(pt));
           dot((v-X(pt)),dXdeta(pt))];

nptsXi = 10;
nptsEta = 10;
npts = nptsXi*nptsEta;

xi_arr = linspace(bnd(1,1),bnd(1,2),nptsXi);
eta_arr = linspace(bnd(2,1),bnd(2,2),nptsEta);
xi_arr = copyVector(xi_arr,nptsEta,1);
eta_arr = copyVector(eta_arr,nptsXi,2);

d_I = norm2((repmat(v,1,npts)-X([xi_arr,eta_arr].')).');
[d,I] = min(d_I);

noItrs = 100;
if d < 1e-13
    return
end
[extrema,nrItr] = newtonsMethodND(f,@(pt)jac(pt,v),[xi_arr(I),eta_arr(I)].',noItrs,1e-15,bnd);

d = norm(v-X(extrema));

function J = jac(pt,x)

J = zeros(2);
J(1,1) = -norm(dXdxi(pt))^2 + dot((x-X(pt)),d2Xdxi2(pt));
J(1,2) = -dot(dXdxi(pt),dXdeta(pt)) + dot((x-X(pt)),d2Xdxideta(pt));
J(1,2) = J(2,1);
J(2,2) = -norm(dXdeta(pt))^2 + dot((x-X(pt)),d2Xdeta2(pt));
end
end