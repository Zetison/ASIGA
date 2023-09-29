function x = abcdeFormula(C)
% Solves polynomial equation up to order 4

x = zeros(size(C,1),size(C,2)-1);
switch size(C,2)
    case 2 % linear equation
        a = C(:,1);
        b = C(:,2);
        indices = a ~= 0;
        x(~indices,1) = NaN;
        
        a = a(indices);
        b = b(indices);
        x(indices) = -b./a;
    case 3 % https://en.wikipedia.org/wiki/Quadratic_formula
        a = C(:,1);
        b = C(:,2);
        c = C(:,3);
        indices = a ~= 0;
        x(~indices,1) = abcdeFormula(C(~indices,2:3));
        x(~indices,2) = NaN;
        
        a = a(indices);
        b = b(indices);
        c = c(indices);
        
        x(indices,1) = (-b + sqrt(b.^2-4*a.*c))./(2*a);
        x(indices,2) = (-b - sqrt(b.^2-4*a.*c))./(2*a);
    case 4 % https://en.wikipedia.org/wiki/Cubic_function
        a = C(:,1);
        b = C(:,2);
        c = C(:,3);
        d = C(:,4);
        
        indices = a ~= 0;
        x(~indices,1:2) = abcdeFormula(C(~indices,2:4));
        x(~indices,3) = NaN;
        
        a = a(indices);
        b = b(indices);
        c = c(indices);
        d = d(indices);
        
        Delta0 = b.^2 - 3*a.*c;
        Delta1 = 2*b.^3-9*a.*b.*c+27*a.^2.*d;
        indices2 = Delta0 ~= 0;
        Ct = zeros(size(C,1),1);
        Ct(indices2) = ((Delta1(indices2) + sqrt(Delta1(indices2).^2 - 4*Delta0(indices2).^3))/2).^(1/3);
        Ct(~indices2) = (Delta1(~indices2)).^(1/3);
        
        xi = -0.5+0.5*sqrt(3)*1i;
        x(indices,1) = -(b+Ct+Delta0./Ct)./(3*a);
        x(indices,2) = -(b+xi*Ct+Delta0./(xi*Ct))./(3*a);
        x(indices,3) = -(b+xi^2*Ct+Delta0./(xi^2*Ct))./(3*a);
    case 5 % http://mathworld.wolfram.com/QuarticEquation.html
        a = C(:,1);
        b = C(:,2);
        c = C(:,3);
        d = C(:,4);
        e = C(:,5);
        
        indices = a ~= 0;
        x(~indices,1:3) = abcdeFormula(C(~indices,2:5));
        x(~indices,4) = NaN;
        
        a = a(indices);
        b = b(indices);
        c = c(indices);
        d = d(indices);
        e = e(indices);
        
        a3 = b./a;
        a2 = c./a;
        a1 = d./a;
        a0 = e./a;

        y = abcdeFormula([ones(size(a2)), -a2, a1.*a3-4*a0, 4*a2.*a0-a1.^2-a3.^2.*a0]);
        y1 = y(:,1);
        R = sqrt(a3.^2/4 - a2 + y1);
        indices2 = R==0;
        D = zeros(size(C,1),1);
        D(~indices2) = sqrt(3/4*a3(~indices2).^2 - R(~indices2).^2 - 2*a2(~indices2) + (4*a3(~indices2).*a2(~indices2)-8*a1(~indices2)-a3(~indices2).^3)./(4*R(~indices2)));
        D(indices2) = sqrt(3/4*a3(indices2).^2 - 2*a2(indices2) + 2*sqrt(y1(indices2).^2-4*a0(indices2)));
        E = zeros(size(C,1),1);
        E(~indices2) = sqrt(3/4*a3(~indices2).^2 - R(~indices2).^2 - 2*a2(~indices2) - (4*a3(~indices2).*a2(~indices2)-8*a1(~indices2)-a3(~indices2).^3)./(4*R(~indices2)));
        E(indices2) = sqrt(3/4*a3(indices2).^2 - 2*a2(indices2) - 2*sqrt(y1(indices2).^2-4*a0(indices2)));
        
        x(indices,1) = -a3/4+R/2+D/2;
        x(indices,2) = -a3/4+R/2-D/2;
        x(indices,3) = -a3/4-R/2+E/2;
        x(indices,4) = -a3/4-R/2-E/2;
        
        
end



