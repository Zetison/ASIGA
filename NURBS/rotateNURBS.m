function nurbs = rotateNURBS(varargin)
options = struct('rotAxis', [0,0,1], ...   % knot vector of revolution parametric direction
                 'x_0', [0,0,0], ...   % knot vector of revolution parametric direction
                 'theta', pi/2);     % angle of revolution in degrees
nurbs = varargin{1};
if nargin > 1
    if numel(varargin) > 2
        newOptions = varargin(2:end);
    else
        newOptions = varargin{2};
    end
    options = updateOptions(options,newOptions);
end
theta = options.theta;
rotAxis = options.rotAxis;
nurbs = ensure3DNURBS(nurbs);
if theta == 0
    return
end
x_0 = options.x_0;
if any(x_0 ~= 0)
    nurbs = translateNURBS(nurbs,-x_0);
end
for patch = 1:numel(nurbs)
    coeffs = nurbs{patch}.coeffs;
    R_x = rotationMatrix(theta, rotAxis);
    
    coeffs = subasgnArr(coeffs,matrixArrayProd(R_x,slc(coeffs, 1:3)),1:3);

    nurbs{patch}.coeffs = coeffs;
end
if any(x_0 ~= 0)
    nurbs = translateNURBS(nurbs,x_0);
end