function nurbs = getPrismShellData(varargin)

options = struct('L', [1,1,1], ...
                 't', 0.1, ...
                 'x_0', [0,0,0], ...
                 't_extraKnots', [], ...
                 'addPMLinfo', false);
                   
if nargin > 0
    if numel(varargin) > 1
        newOptions = varargin;
    else
        newOptions = varargin{1};
    end
    options = updateOptions(options,newOptions);
end

% If the centerCoordinate is not given, the prism is placed with its center
% at the origin

t = options.t;
t_e = options.t_extraKnots;
L = options.L;
Ltot = L+2*t;
nurbs = getPrismData('x_0',options.x_0,'L',Ltot);
nurbs = insertKnotsInNURBS(nurbs,{[t,t+L(1)]/Ltot(1),[t,t+L(2)]/Ltot(2),[t,t+L(3)]/Ltot(3)});
nurbs = explodeNURBS(nurbs,1:3,1);
for i = 1:numel(t_e)
    nurbs(2:3:end) = insertKnotsInNURBS(nurbs(2:3:end),{[t_e(i),L(1)-t_e(i)]/L(1),[],[]});
    nurbs([4:6,13:15,22:24]) = insertKnotsInNURBS(nurbs([4:6,13:15,22:24]),{[],[t_e(i),L(2)-t_e(i)]/L(2),[]});
    nurbs(10:18) = insertKnotsInNURBS(nurbs(10:18),{[],[],[t_e(i),L(3)-t_e(i)]/L(3)});
end

nurbs(1:9) = flipNURBSparametrization(nurbs(1:9),3);
nurbs([1:3,10:12,19:21]) = flipNURBSparametrization(nurbs([1:3,10:12,19:21]),2);
nurbs(1:3:end) = flipNURBSparametrization(nurbs(1:3:end),1);
nurbs(1) = permuteNURBS(nurbs(1),[1,3,2]);
nurbs([5:6,8:9,11:13,16,20:22,25]) = permuteNURBS(nurbs([5:6,8:9,11:13,16,20:22,25]),[2,1,3]);

if options.addPMLinfo
    for i = 1:numel(nurbs)
        switch i
            case {1,3,7,9,19,21,25,27}
                nurbs{i}.isPML = [1,1,1];  
            case {10,12,16,18}
                nurbs{i}.isPML = [1,1,0];  
            case {2,6,22,26}
                nurbs{i}.isPML = [0,1,1]; 
            case {4,8,20,24}
                nurbs{i}.isPML = [1,0,1];   
            case {5,23}
                nurbs{i}.isPML = [0,0,1];    
            case {13,17}
                nurbs{i}.isPML = [0,1,0];   
            case {11,15}
                nurbs{i}.isPML = [1,0,0];    
        end
    end
end
% figure,plotNURBS(nurbs,'resolution',[10,10,10],'plotParmDir',true)
nurbs(14) = [];