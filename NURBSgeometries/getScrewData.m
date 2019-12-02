function nurbs = getScrewData(r,R,l,h,t,type,s)


Eta = [0 0 1 1];
Zeta = [0 0 0.5 1 1];
if type == 1
    Xi = [0 0 0 0.5 0.5 1 1 1];
    controlPts = zeros(4,5,2,3);
    theta = pi/4;
    x = r*cos(theta);
    y = 2*x - r;
    w = cos(theta/2);
    controlPts(:,1,1,1) = [r, 0, 0, 1];
    controlPts(:,2,1,1) = [r, y, 0, w];
    controlPts(:,3,1,1) = [x, x, 0, 1];
    controlPts(:,4,1,1) = [y, r, 0, w];
    controlPts(:,5,1,1) = [0, r, 0, 1];
    
    controlPts(:,1,2,1) = [l, 0, 0, 1];
    controlPts(:,2,2,1) = [l, l/2, 0, 1];
    controlPts(:,3,2,1) = [l, l, 0, 1];
    controlPts(:,4,2,1) = [l/2, l, 0, 1];
    controlPts(:,5,2,1) = [0, l, 0, 1];
    
    controlPts(:,1,1,2) = [r, 0, h, 1];
    controlPts(:,2,1,2) = [r, y, h, w];
    controlPts(:,3,1,2) = [x, x, h, 1];
    controlPts(:,4,1,2) = [y, r, h, w];
    controlPts(:,5,1,2) = [0, r, h, 1];
    
    controlPts(:,1,2,2) = [l, 0, h, 1];
    controlPts(:,2,2,2) = [l, l/2, h, 1];
    controlPts(:,3,2,2) = [l, l, h, 1];
    controlPts(:,4,2,2) = [l/2, l, h, 1];
    controlPts(:,5,2,2) = [0, l, h, 1];
    
    x = R*cos(theta);
    y = 2*x - R;
    
    controlPts(:,1,1,3) = [R, 0, t, 1];
    controlPts(:,2,1,3) = [R, y, t, w];
    controlPts(:,3,1,3) = [x, x, t, 1];
    controlPts(:,4,1,3) = [y, R, t, w];
    controlPts(:,5,1,3) = [0, R, t, 1];
    
    controlPts(:,1,2,3) = [l, 0, t, 1];
    controlPts(:,2,2,3) = [l, l/2, t, 1];
    controlPts(:,3,2,3) = [l, l, t, 1];
    controlPts(:,4,2,3) = [l/2, l, t, 1];
    controlPts(:,5,2,3) = [0, l, t, 1];
    
else
    Xi = [0 0 0 1/3 1/3 2/3 2/3 1 1 1];
    controlPts = zeros(4,7,2,3);
    theta = pi/6;
    x = r*cos(theta);
    y = r*sin(theta);
    w = cos(theta/2);
    y1 = y - cot(theta)*(r-x);
    % Y = y - cot(theta)*(X-x)
    % Y = x - cot(2*theta)*(X-y)
    % x - cot(2*theta)*(X-y) = y - cot(theta)*(X-x) 
    % => (1-cot(theta))*x - (1-cot(2*theta))*y = (cot(2*theta)-cot(theta))*X 
    x1 = ((1-cot(theta))*x - (1-cot(2*theta))*y)/(cot(2*theta)-cot(theta));
    controlPts(:,1,1,1) = [r, 0, 0, 1];
    controlPts(:,2,1,1) = [r, y1, 0, w];
    controlPts(:,3,1,1) = [x, y, 0, 1];
    controlPts(:,4,1,1) = [x1, x1, 0, w];
    controlPts(:,5,1,1) = [y, x, 0, 1];
    controlPts(:,6,1,1) = [y1,r, 0, w];
    controlPts(:,7,1,1) = [0, r, 0, 1];
    
    controlPts(:,1,2,1) = [l, 0, 0, 1];
    controlPts(:,2,2,1) = [l, s/2, 0, 1];
    controlPts(:,3,2,1) = [l, s, 0, 1];
    controlPts(:,4,2,1) = [(s+l)/2, (s+l)/2, 0, 1];
    controlPts(:,5,2,1) = [s, l, 0, 1];
    controlPts(:,6,2,1) = [s/2, l, 0, 1];
    controlPts(:,7,2,1) = [0, l, 0, 1];
    
    controlPts(:,1,1,2) = [r, 0, h, 1];
    controlPts(:,2,1,2) = [r, y1,h, w];
    controlPts(:,3,1,2) = [x, y, h, 1];
    controlPts(:,4,1,2) = [x1, x1, h, w];
    controlPts(:,5,1,2) = [y, x, h, 1];
    controlPts(:,6,1,2) = [y1,r, h, w];
    controlPts(:,7,1,2) = [0, r, h, 1];
    
    controlPts(:,1,2,2) = [l, 0, h, 1];
    controlPts(:,2,2,2) = [l, s/2, h, 1];
    controlPts(:,3,2,2) = [l, s, h, 1];
    controlPts(:,4,2,2) = [(s+l)/2, (s+l)/2, h, 1];
    controlPts(:,5,2,2) = [s, l, h, 1];
    controlPts(:,6,2,2) = [s/2, l, h, 1];
    controlPts(:,7,2,2) = [0, l, h, 1];
        
    x = R*cos(theta);
    y = R*sin(theta);
    w = cos(theta/2);
    y1 = y - cot(theta)*(R-x);
    x1 = ((1-cot(theta))*x - (1-cot(2*theta))*y)/(cot(2*theta)-cot(theta));
    
    
    controlPts(:,1,1,3) = [R, 0, t, 1];
    controlPts(:,2,1,3) = [R, y1,t, w];
    controlPts(:,3,1,3) = [x, y, t, 1];
    controlPts(:,4,1,3) = [x1, x1, t, w];
    controlPts(:,5,1,3) = [y, x, t, 1];
    controlPts(:,6,1,3) = [y1,R, t, w];
    controlPts(:,7,1,3) = [0, R, t, 1];
    
    controlPts(:,1,2,3) = [l, 0, t, 1];
    controlPts(:,2,2,3) = [l, s/2, t, 1];
    controlPts(:,3,2,3) = [l, s, t, 1];
    controlPts(:,4,2,3) = [(s+l)/2, (s+l)/2, t, 1];
    controlPts(:,5,2,3) = [s, l, t, 1];
    controlPts(:,6,2,3) = [s/2, l, t, 1];
    controlPts(:,7,2,3) = [0, l, t, 1];
end

% if nargin == 4
%     for i = 1:2
%         for j = 1:2
%             for k = 1:2
%                 controlPts(1:3,i,j,k) = controlPts(1:3,i,j,k) + centerCoordinate;
%             end
%         end
%     end
% end

nurbs = createNURBSobject(controlPts,{Xi, Eta, Zeta});
