function nurbs = getEllipsoidData(varargin)
% set default values
options = struct('Xi', [0,0,0,1,1,2,2,3,3,4,4,4]/4, ...  % knot vector for azimuthal direction
                 'Eta', [0,0,0,1,1,2,2,2]/2, ...         % knot vector for polar direction
                 'C', 1, ...
                 'alignWithAxis', 'Zaxis', ...
                 'x_0',[0, 0, 0],...
                 'alpha', 0,...                          % rotation angle about alignWithAxis
                 'theta', 2*pi,...
                 'theta_eta', pi,...
                 'parm', 1, ...
                 't', 0, ...
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

alignWithAxis = options.alignWithAxis;
C = options.C;
t = options.t;
prec = options.prec;
if numel(C) == 1
    C = C*ones(1,3,prec);
elseif numel(C) == 2
    C = [C,C(2)];
end

switch options.parm
    case 1
        optionsArc.Xi = options.Eta;
        optionsArc.theta = options.theta_eta;
        nurbs = getArcData(optionsArc);
        options.theta = 2*pi;
        options.rotAxis = [1,0,0];
        nurbs = rotateNURBS(nurbs,'rotAxis','Yaxis','theta',pi/2);
        nurbs = rotateNURBS(nurbs,'rotAxis','Zaxis','theta',-pi/2);
        nurbs = revolveNURBS(nurbs,'theta',2*pi,'rotAxis','Zaxis','Xi',options.Xi);
        nurbs = permuteNURBS(nurbs,[2,1]);
    case 2
        Xi = [0 0 0 0 0 1 1 1 1 1];
        Eta = [0 0 0 0 0 1 1 1 1 1];

        if strcmp(prec,'sym')
            sr2 = sqrt(sym('2'));                                                        
            sr3 = sqrt(sym('3'));                                                        
            sr6 = sqrt(sym('6'));   
        else
            sr2 = sqrt(2);                                                        
            sr3 = sqrt(3);                                                        
            sr6 = sqrt(6);   
        end
        controlPts = zeros(4,5,5,prec); 
        controlPts(:,1,1,1) = [4*(1-sr3)    4*(1-sr3)       4*(sr3-1)       4*(3-sr3)   ];
        controlPts(:,2,1,1) = [-sr2         sr2*(sr3-4)     sr2*(4-sr3)     sr2*(3*sr3-2)];
        controlPts(:,3,1,1) = [0            4*(1-2*sr3)/3   4*(2*sr3-1)/3	4*(5-sr3)/3]; 
        controlPts(:,2,2,1) = [-(3*sr3-2)/2	(2-3*sr3)/2     (sr3+6)/2       (sr3+6)/2];
        controlPts(:,3,2,1) = [0            sr2*(2*sr3-7)/3	5*sr6/3         sr2*(sr3+6)/3];
        controlPts(:,3,3,1) = [0            0               4*(5-sr3)/3     4*(5*sr3-1)/9];
        pairs = [1,2; 1,3; 2,3];
        for l = 1:3
            i = pairs(l,1);
            j = pairs(l,2);
            controlPts(1,i,j) = controlPts(2,j,i);
            controlPts(2,i,j) = controlPts(1,j,i);
            controlPts(3:4,i,j) = controlPts(3:4,j,i);
        end
        for i = 4:5
            for j = 1:3
                controlPts(1,i,j) = -controlPts(1,6-i,j);
                controlPts(2:4,i,j) = controlPts(2:4,6-i,j);
            end
        end
        for i = 1:5
            for j = 4:5
                controlPts(2,i,j) = -controlPts(2,i,6-j);
                controlPts([1,3,4],i,j) = controlPts([1,3,4],i,6-j);
            end
        end
        maxWeight = 4*(3-sr3);
        for i = 1:5
            for j = 1:5
                controlPts(1:3,i,j) = controlPts(1:3,i,j)/controlPts(4,i,j);
                controlPts(4,i,j) = controlPts(4,i,j)/maxWeight;
            end
        end  
        nurbs = cell(1,6);
        nurbs(1) = createNURBSobject(controlPts,{Xi, Eta});

        for i = 2:4
            nurbs(i) = rotateNURBS(nurbs(1),'theta',(i-1)*pi/2,'rotAxis','Yaxis');
        end
        nurbs(5) = rotateNURBS(nurbs(1),'theta',pi/2,'rotAxis','Xaxis');
        nurbs(6) = rotateNURBS(nurbs(1),'theta',-pi/2,'rotAxis','Xaxis');
end
    
switch alignWithAxis
    case 'Xaxis' 
        nurbs = rotateNURBS(nurbs,'rotAxis',[1,1,1],'theta',120*pi/180);
    case 'Yaxis'
        nurbs = rotateNURBS(nurbs,'rotAxis',[1,1,1],'theta',2*120*pi/180);
end
nurbs = rotateNURBS(nurbs,'theta',options.alpha,'rotAxis',alignWithAxis);
if t > 0
    nurbs_i = scaleNURBS(nurbs,C-t);
    nurbs_o = scaleNURBS(nurbs,C);
    nurbs = loftNURBS({nurbs_i,nurbs_o});
else
    nurbs = scaleNURBS(nurbs,C);
end
nurbs = translateNURBS(nurbs,options.x_0);
if options.uniformDegree
    nurbs = makeUniformNURBSDegree(nurbs);
end
