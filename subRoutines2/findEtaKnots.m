newEtaValues1 = zeros(1,mesh1);
theta_prev = theta_eta1;

for i = 1:mesh1
    fun = @(theta) findArcLength(R_a,Upsilon,theta_prev,theta) - totArcLength1/(mesh1+1);
    theta = bisection(fun, theta_prev, pi, 100, 1e-12);
    
    phi = 0;
    x_vec = [sqrt(R_a^2-Upsilon^2)*sin(theta)*cos(phi);
             sqrt(R_a^2-Upsilon^2)*sin(theta)*sin(phi);
             R_a*cos(theta)];
    fun = @(eta) norm(x_vec - evaluateNURBS(fluid, [0, eta, 1]));
    
    newEtaValues1(i) = fminsearchbnd(fun,eta1/2,0,eta1);
    theta_prev = theta;
end

theta_prev = 0;
newEtaValues1 = fliplr(newEtaValues1);

newEtaValues2 = zeros(1,mesh3);
for i = 1:mesh3
    fun = @(theta) findArcLength(R_a,Upsilon,theta_prev,theta) - totArcLength2/(mesh3+1);
    theta = bisection(fun, theta_prev, theta_eta2, 100, 1e-12);
    
    phi = 0;
    x_vec = [sqrt(R_a^2-Upsilon^2)*sin(theta)*cos(phi);
             sqrt(R_a^2-Upsilon^2)*sin(theta)*sin(phi);
             R_a*cos(theta)];
    fun = @(eta) norm(x_vec - evaluateNURBS(fluid, [0, eta, 1]));
    
    newEtaValues2(i) = fminsearchbnd(fun,(1+eta2)/2,eta2,1);
    theta_prev = theta;
end
newEtaValues2 = fliplr(newEtaValues2);