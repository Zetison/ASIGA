function nurbs = getArcData(varargin)
options = struct('Xi', [0,0,0,1,1,2,2,3,3,3]/3, ...   % knot vector of revolution parametric direction
                 'R', 1, ...   % Radius
                 'theta', 2*pi);     % angle of revolution in degrees
if nargin > 0
    if numel(varargin) > 1
        newOptions = varargin;
    else
        newOptions = varargin{1};
    end
    options = updateOptions(options,newOptions);
end
Xi = options.Xi;
theta = options.theta;
if Xi(end) ~= 1
    Xi = Xi/Xi(end);
end
n_xi = numel(Xi)-3;
controlPts = zeros(3,n_xi);
controlPts(:,1)   = [1, 0, 1];

counter = 1;
theta_prev = 0;

i = 4;
while i < length(Xi)
    P1 = controlPts(1:2,counter);
    
    i_prev = i;
    while sum(Xi(i) == Xi) < 2
        i = i + 1;
    end
    theta_i = Xi(i)*theta;
    P3 = [cos(theta_i); sin(theta_i)];
    t1 = (1-P1(1)*P3(1)-P1(2)*P3(2))/(P1(2)*P3(1)-P1(1)*P3(2));
    P2 = P1 + t1*[P1(2); -P1(1)];
    w2 = cos((theta_i-theta_prev)/2);
    
    controlPtsTemp = zeros(3,3);
    controlPtsTemp(:,1) = [P1; 1];
    controlPtsTemp(:,2) = [P2; w2];
    controlPtsTemp(:,3) = [P3; 1];
    Xi_temp = [0,0,0,1,1,1];
    nurbs_temp = createNURBSobject(controlPtsTemp,{Xi_temp});
    
    nurbs_temp = insertKnotsInNURBS(nurbs_temp,{linspace2(0,1,i-i_prev)});
    controlPts(:,counter+1:counter+2+i-i_prev) = nurbs_temp{1}.coeffs(:,2:end);
    counter = counter + 2+i-i_prev;

    theta_prev = theta_i;
    i = i + 2;
end
nurbs = createNURBSobject(controlPts,Xi);
nurbs = scaleNURBS(nurbs,options.R);
