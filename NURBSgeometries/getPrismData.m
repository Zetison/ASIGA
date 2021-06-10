function nurbs = getPrismData(varargin)

options = struct('X', NaN, ...
                 'L', [1,1,1], ...
                 'x_0', [0,0,0]);
                   
if nargin > 0
    if numel(varargin) > 1
        newOptions = varargin;
    else
        newOptions = varargin{1};
    end
    options = updateOptions(options,newOptions);
end

% If the centerCoordinate is not given, the prism is placed with its center
% at the options.L/2

Xi = [0 0 1 1];
Eta = [0 0 1 1];
Zeta = [0 0 1 1];

coeffs = ones(4,2,2,2);
X = options.X;
if isnan(X)
    L = options.L;
    if numel(L) == 1
        L = L*ones(1,3);
    end
    X = [  0,   0,   0;
           L(1),0,   0;
           0,   L(2),0;
           L(1),L(2),0;
           0,   0,   L(3);
           L(1),0,   L(3);
           0,   L(2),L(3);
           L(1),L(2),L(3)];
end
coeffs(1:3,:,:,:) = reshape(X.',3,2,2,2);

nurbs = createNURBSobject(coeffs,{Xi, Eta, Zeta});
nurbs = translateNURBS(nurbs,options.x_0);
