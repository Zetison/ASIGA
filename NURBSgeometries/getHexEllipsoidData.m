function nurbs = getHexEllipsoidData(varargin)
% set default values
options = struct('C', 1, ...
                 'alignWithAxis', 'Zaxis', ...
                 'x_0',[0, 0, 0],...
                 'prec', 'double', ...
                 'uniformDegree', true);
             
if nargin > 0
    if numel(varargin) > 1
        newOptions = varargin;
    else
        newOptions = varargin{1};
    end
    options = updateOptions(options,newOptions);
end
options.t = 0;
options.parm = 2;
nurbsSurf = getEllipsoidData(options);
nurbs = cell(1,6);
nurbs{1} = orientNURBS(nurbsSurf{4},4);
nurbs{2} = orientNURBS(nurbsSurf{2},5);
nurbs{3} = orientNURBS(nurbsSurf{5},4);
nurbs{4} = orientNURBS(nurbsSurf{6},6);
nurbs{5} = orientNURBS(nurbsSurf{3},1);
nurbs{6} = orientNURBS(nurbsSurf{1},0);
nurbs = GordonHall(nurbs);