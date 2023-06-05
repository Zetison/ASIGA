function nurbs = getRectangleData(varargin)

options = struct('X', NaN, ...
                 'L', [1,1], ...
                 'x_0', [0,0]);
                   
if nargin > 0
    if numel(varargin) > 1
        newOptions = varargin;
    else
        newOptions = varargin{1};
    end
    options = updateOptions(options,newOptions);
end

% If the x_0 (center coordinate) is not given, the rectangle is placed with its center
% at the options.L/2

Xi = [0 0 1 1];
Eta = [0 0 1 1];

coeffs = ones(3,2,2);
X = options.X;
if isnan(X)
    L = options.L;
    if numel(L) == 1
        L = L*ones(1,2);
    end
    X = [  0,   0;
           L(1),0;
           0,   L(2);
           L(1),L(2)];
end
coeffs(1:2,:,:) = reshape(X.',2,2,2);

nurbs = createNURBSobject(coeffs,{Xi, Eta});
nurbs = translateNURBS(nurbs,options.x_0);
