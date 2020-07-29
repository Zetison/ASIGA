function nurbs = revolveNURBS(varargin)
% set default values
options = struct('Xi', [0,0,0,1,1,2,2,3,3,3]/3, ...   % knot vector of revolution parametric direction
                 'theta', 2*pi, ...             % angle of revolution in degrees
                 'rotAxis', [0, 0, 1],...       % rotation axis
                 'x_1', [0, 0, 0]);             % origin of rotation axis
nurbs = varargin{1};
if nargin > 1
    if numel(varargin) > 2
        newOptions = varargin(2:end);
    else
        newOptions = varargin{2};
    end
    options = updateOptions(options,newOptions);
end
nurbs = ensure3DNURBS(nurbs);
if isa(options.rotAxis,'char')
    switch options.rotAxis
        case 'Xaxis'
            options.rotAxis = [1, 0, 0];
        case 'Yaxis'
            options.rotAxis = [0, 1, 0];
        case 'Zaxis'
            options.rotAxis = [0, 0, 1];
    end
end
x_1 = options.x_1;
Xi = options.Xi;
n = options.rotAxis;
n = n/norm(n);  % make sure n is a unit vector
nurbsArc = ensure3DNURBS(getArcData(strctCp(options, {'Xi','theta'})));
controlPts = nurbsArc{1}.coeffs(1:3,:);
wArc = nurbsArc{1}.coeffs(end,:);
[~, A] = orthogonalTransform(NaN(1,3), n);
controlPts = (A*controlPts).';
n_xi = numel(Xi)-(2+1);
for i = 1:numel(nurbs)
    sizes = size(nurbs{i}.coeffs);
    number = nurbs{i}.number;
    d = nurbs{i}.d;
    x = nurbs{i}.coeffs(1:d,:);
    w = nurbs{i}.coeffs(d+1,:);
    x = reshape(x,d,prod(sizes(2:end))).';
    t = (x-x_1)*n.';
    r = x_1 + t*n; % parametrization of rotation axis
    R = norm2(r-x); % distance from parametrization of rotation axis to control point
    coeffs = zeros([sizes,n_xi]);
    for j = 1:n_xi
        x = r + R*controlPts(j,:);
        
        subs = {[{1:d},repmat({':'},1,ndims(coeffs)-2),{j}]};
        coeffs = subasgnArr(coeffs,reshape(x.',[d,number]),subs);
        subs{1}{1} = d+1;
        coeffs = subasgnArr(coeffs,reshape(wArc(j)*w,[1,number]),subs);
    end
    nurbs{i}.coeffs = coeffs;
    nurbs{i}.d_p = nurbs{i}.d_p + 1;
    nurbs{i}.number = [sizes(2:end),n_xi];
    nurbs{i}.degree = [nurbs{i}.degree,2];
    nurbs{i}.knots = [nurbs{i}.knots(:)', {Xi}];
end
