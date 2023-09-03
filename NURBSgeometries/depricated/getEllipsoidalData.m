function nurbs = getEllipsoidalData(newOptions)
error('Depricated: use getEllipsoidData() instead')
% 
% %% Interpret input arguments
% % set default values
% options = struct('R', [1,1,1], ...
%                  'alignWithAxis', 'Xaxis', ...
%                  'x_0',[0, 0, 0],...
%                  'alpha', 0,...
%                  'parm', 1, ...
%                  't', 0, ...
%                  'iXi', [1,1,2,2,3,3]/4,...
%                  'iEta', [1,1]/2);
%              
% options = updateOptions(options,newOptions);
% 
% alignWithAxis = options.alignWithAxis;
% iXi = options.iXi;
% iEta = options.iEta;
% R = options.R;
% t = options.t;
% if numel(R) == 1
%     R = R*ones(1,3);
% elseif numel(R) == 2
%     R = [R,R(2)];
% end
% nurbs = getUnitSphereData(options.parm,class(R),iXi,iEta);
% 
%     
% switch alignWithAxis
%     case 'Yaxis' 
%         nurbs = rotateNURBS(nurbs,120*pi/180,[1,1,1]);
%     case 'Zaxis'
%         nurbs = rotateNURBS(nurbs,2*120*pi/180,[1,1,1]);
% end
% nurbs = rotateNURBS(nurbs,options.alpha,alignWithAxis);
% if t > 0
%     nurbs_i = scaleNURBS(nurbs,R-t);
%     nurbs_o = scaleNURBS(nurbs,R);
%     nurbs = linearSurfaceToVolume(nurbs_o,nurbs_i);
% else
%     nurbs = scaleNURBS(nurbs,R);
% end
% nurbs = translateNURBS(nurbs,options.x_0);
% makeUniformNURBSDegree