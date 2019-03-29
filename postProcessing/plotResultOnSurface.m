close all
Xi = varCol.nurbs.knots{1};
Eta = varCol.nurbs.knots{2};
p_xi = varCol.nurbs.degree(1);
p_eta = varCol.nurbs.degree(2);

weights = varCol.weights;
controlPts = varCol.controlPts;
noDofs = varCol.noDofs;
nurbs = varCol.nurbs;
d = varCol.dimension;
runInParallell = varCol.runInParallell;
omega = varCol.omega;
model = varCol.model;

xi = 0.25;
N = 1000;
eta_arr = linspace(0,0.5,N);
u_arr = zeros(size(eta_arr));
x_arr = zeros(size(eta_arr));
for i = 1:N
    eta = eta_arr(i);
    u_arr(i) = numericalSolEval_final_surf(xi, eta, p_xi, p_eta, Xi, Eta, weights, controlPts, U_fluid_o);
    v = evaluateNURBS(fluid,[xi,eta]);
    x_arr(i) = v(1);
end

plot(eta_arr,real(u_arr))

title(sprintf('%s: real part of the solution along xi = const',saveName),'interpreter','none')

ylabel('real(u)')
xlabel('eta - parameter')
% save('test.mat', 'x_arr', 'u_arr')