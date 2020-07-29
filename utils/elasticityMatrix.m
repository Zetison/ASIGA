function C = elasticityMatrix(E,nu,d)
if nargin < 3
    d = 3;
end
switch d
    case 3
        C = zeros(6,6);
        C(1:3,1:3) = E/(1+nu)/(1-2*nu)*...
                     [1-nu nu nu;
                      nu 1-nu nu; 
                      nu nu 1-nu];
        C(4:6,4:6) = 0.5*E/(1+nu)*eye(3);
    case 2
    	C = 1/(1-nu^2)*[1, nu, 0;
                        nu, 1, 0;
                        0,  0, (1-nu)/2]; % E0 is factored out from C
end