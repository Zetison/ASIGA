function [nurbs, nurbs2, nurbs3, nurbs4, nurbs5] = getBCDataParts(a, b, l, g2, g3, alpha, beta, s, c)
h = g2*tan(alpha/2); % = g2/tan((pi-alpha)/2)

x1 = tan(alpha/2)*(g2/tan(alpha) + h);
x2 = g3*tan(alpha);

lowerCpts = [  -l-g2-g3,    0,  0,          1; % 0
               -l-g2-g3,    0, (-b+h+x2)/2, 1;
               -l-g2-g3,    0, -b+h+x2,     1; % 1
               -l-g2-g3/2,  0, -b+h+x2/2,   1;
               -l-g2,       0, -b+h,        1; % 2
               -l-x1,       0, -b,          cos(alpha/2);
               -l,          0, -b,          1; % 3
               -l/2,        0, -b,          1;
                0,          0, -b,          1; % 4
                a,          0, -b,          1/sqrt(2);
                a,          0,  0,          1];% 5
            
righ_sCpts1 = lowerCpts;
righ_sCpts1(:,2) = abs(lowerCpts(:,3))/tan(pi/2-beta/4);
righ_sCpts1(:,4) = righ_sCpts1(:,4)*cos(beta/4);


R_x = rotationXaxis(beta/2);

righ_sCpts2 = lowerCpts;
righ_sCpts2(:,1:3) = (R_x*righ_sCpts2(:,1:3).').';

coeffs = zeros(4,3,11);
coeffs(:,1,:) = lowerCpts.';
coeffs(:,2,:) = righ_sCpts1.';
coeffs(:,3,:) = righ_sCpts2.';

Xi = [0 0 0 1 1 2 2 2]/2;
Eta = [0 0 0 1 1 2 2 3 3 4 4 5 5 5]/5;

nurbs_lower = createNURBSobject(coeffs,{Xi, Eta});

coeffs = zeros(4,5,11);
coeffs(:,3:5,:) = nurbs_lower.coeffs;
coeffs(:,1:2,:) = nurbs_lower.coeffs(:,3:-1:2,:);
coeffs(2,1:2,:) = -coeffs(2,1:2,:);
coeffs2 = coeffs(:,:,5:end);
nurbs = createNURBSobject(coeffs2,{Xi, [0 0 0 1 1 2 2 3 3 3]/3});

coeffs = coeffs(:,:,1:5);
coeffs(:,7,:) = coeffs(:,1,:);
coeffs(:,6,:) = coeffs(:,2,:);
for i = 1:5
    coeffs(3,6,i) = norm(coeffs(2:3,6,i));
end
coeffs(2,6,:) = 0;
Xi = [0 0 0 1 1 2 2 3 3 3]/3;
nurbs2 = createNURBSobject(coeffs,{Xi, [0,0,0,1,1,2,2,2]/2});


%% Upper main body
C_4 = c + b*cos(beta/2);
C_3 = b*sin(beta/2)-s;
C_2 = (2*C_4+C_3*tan(beta/2))/C_3^3;
C_1 = -(3*C_4+C_3*tan(beta/2))/C_3^2;

P = @(y) c + C_1*(y-s).^2 + C_2*(y-s).^3;

npts = 6;
dy = C_3/(npts-1);

coeffs = zeros(4,11,3);

y = s + (npts-1)*dy;
for i = 1:npts
    line = zeros(3,4);
    line(1,:)   = [-l,   y, P(y), 1];
    line(2,:)   = [-l/2, y, P(y), 1];
    line(3,:) = [0,    y, P(y), 1];
    coeffs(:,2*i-1,:) = line.';
    y = y - dy;
end

t1 = 120*pi/180;
p2 = t1/12;
w2 = cos(p2/2);
for i = 1:5
    coeffs(1:3,2*i,:) = (coeffs(1:3,2*i-1,:)+coeffs(1:3,2*i+1,:))/2;
    coeffs(4,2*i,[1,2,3]) = w2;
end
nurbs3 = createNURBSobject(coeffs,{[0,0,0,1,1,2,2,3,3,4,4,5,5,5]/5, [0,0,0,1,1,1]});
coeffs2 = coeffs;
coeffs2(2,:,:) = -coeffs(2,:,:);
nurbs4 = createNURBSobject(coeffs2,{[0,0,0,1,1,2,2,3,3,4,4,5,5,5]/5, [0,0,0,1,1,1]});

coeffs3 = zeros(4,2,3);
coeffs3(:, 1, :) = coeffs(:, end, :);
coeffs3(:, 2, :) = coeffs2(:, end, :);
nurbs5 = createNURBSobject(coeffs3,{[0,0,1,1], [0,0,0,1,1,1]});


function R_x = rotationXaxis(alpha)

R_x = [1, 0,           0;
       0, cos(alpha), -sin(alpha);
       0, sin(alpha),  cos(alpha)];
