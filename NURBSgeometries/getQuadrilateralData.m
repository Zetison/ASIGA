function nurbs = getQuadrilateralData(varargin)
options = struct('X', [0,0,0;
                       1,0,0;
                       0,1,0;
                       1,1,0]);   % Corners of quadrilateral in 3D
if nargin > 0
    if numel(varargin) > 1
        newOptions = varargin;
    else
        newOptions = varargin{1};
    end
    options = updateOptions(options,newOptions);
end
X = options.X.';
d = size(X,1);
Xi = [0 0 1 1];
Eta = [0 0 1 1];
coeffs = ones(d+1,2,2);
coeffs(1:d,:,:) = reshape(X,d,2,2);

nurbs = createNURBSobject(coeffs,{Xi, Eta});
