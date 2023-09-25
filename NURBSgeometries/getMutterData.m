function nurbs = getMutterData(varargin)
options = struct('R1', 6,...
                 'R2', 6.5,...
                 'R3', 9,...
                 't1', 9.8, ...
                 't2', 10.8, ...
                 't3', 10.8-2*0.50675752076);
if nargin > 0
    if numel(varargin) > 1
        newOptions = varargin;
    else
        newOptions = varargin{1};
    end
    options = updateOptions(options,newOptions);
end
R1 = options.R1;
R2 = options.R2;
R3 = options.R3;
t1 = options.t1;
t2 = options.t2;
t3 = options.t3;
options.d_p = 2;
options.X = [0,0,-R3;
             0,t2,-R3;
             0,0,-R2;
             0,t2,-R2];
nurbs1 = getPrismData(options);
delta = (t2-t1)/2;
options.X = [options.X(3:4,:);
             0,delta,-R1;
             0,t2-delta,-R1];
nurbs1(2) = getPrismData(options);
twelfth = 2*pi/12;
nurbs1 = revolveNURBS(nurbs1,'Xi', [0,0,0,1,1,1], 'theta', -twelfth, 'rotAxis', [0, 1, 0]);
nurbs1 = elevateNURBSdegree(nurbs1,[0,0,1]);
nurbs1 = insertKnotsInNURBS(nurbs1,{[],[],0.5});
delta2 = (t2-t3)/2;
coeffs = [0,0,-R3, 1;
          R3*tan(twelfth),delta2,-R3, 1].';
nurbs_bottomCurve = createNURBSobject(coeffs,[0,0,1,1]);
nurbs_bottomCurve = elevateNURBSdegree(nurbs_bottomCurve,2);
nurbs_bottomCurve = insertKnotsInNURBS(nurbs_bottomCurve,{0.5});
nurbs_bottomCurve{1}.coeffs(4,:) = nurbs1{1}.coeffs(4,1,1,:);
nurbs_bottomCurve{1}.coeffs(1:3,2) = nurbs1{1}.coeffs(1:3,1,1,2);
nurbs_bottomCurve{1}.coeffs(2,3) = delta2*0.094458282367901/0.50675752076;
nurbs_bottomCurve{1}.coeffs(2,4) = delta2*0.350871810563701/0.50675752076;

nurbs_bottomCurve_mirror = translateNURBS(mirrorNURBS(nurbs_bottomCurve,'y'),[0,t2,0]);
nurbs_bottom = loftNURBS({nurbs_bottomCurve,nurbs_bottomCurve_mirror});
nurbs1(3) = loftNURBS({nurbs_bottom,subNURBS(nurbs1(1),'at',[0,0;1,0;0,0])});
nurbs1(3) = permuteNURBS(nurbs1(3),[2,3,1]);

% Ensure smoothness
P1 = nurbs1{3}.coeffs(1:3,1,end,end);
P2 = nurbs1{3}.coeffs(1:3,1,end,end-1);
P3 = nurbs1{3}.coeffs(1:3,1,1,end);
n = cross(P2-P1,P3-P1);
x = nurbs1{3}.coeffs(1:3,1,1,end-1);
nurbs1{3}.coeffs(1,1,1,end-1) = (dot(P1,n) - dot(x(2:3),n(2:3)))/n(1); % Ensure point to lie in the plane spanned by P1, P2 and P3

nurbs1 = elevateNURBSdegree(nurbs1,[0,2,0]);
nurbs1 = permuteNURBS(nurbs1,[2,3,1]);
nurbs2 = mirrorNURBS(nurbs1,'x');
nurbs2 = flipNURBSparametrization(nurbs2,2);
nurbs3 = uniteNURBS({nurbs1,nurbs2});
nurbs = cell(1,36);
for i = 0:5
    nurbs(1+i*6:(i+1)*6) = rotateNURBS(nurbs3,'rotAxis', [0,1,0], 'theta', 2*twelfth*i);
end