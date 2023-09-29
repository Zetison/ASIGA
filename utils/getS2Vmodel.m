function app = getS2Vmodel(app,model)
app.AlgorithmForCase13.Value = 'A13_21';
app.AlgorithmForCase135.Value = 'A135_21';
app.EnforcesymmetryaboutthexyplaneCheckBox.Value = false;
app.EnforcesymmetryaboutthexzplaneCheckBox.Value = false;
app.EnforcesymmetryabouttheyzplaneCheckBox.Value = false;
app.xminEditField.Value = -Inf;
app.xmaxEditField.Value = Inf;
app.yminEditField.Value = -Inf;
app.ymaxEditField.Value = Inf;
app.zminEditField.Value = -Inf;
app.zmaxEditField.Value = Inf;
app.zmaxEditField.Value = Inf;
app.Seth_maxEditField.Value = 0.3;
switch model
    case 'S1'
        app.EpsEditField.Value = 1e-10;
        app.ThicknessEditField.Value = 0.1;
        app.noExtraLayersEditField.Value = 0;
        app.ColorfunctionEditField.Value = '@(x) log10(abs(norm2(x(:,1:3))-1))';
        app.nurbsObjects{1,1}.nurbs = read_g2(['NURBSgeometries/g2files/' model '.g2']);
        app.AnglethresholdEditField.Value = 140; % Threshold for a "sharp" angle
        app.AlgorithmForCase13.Value = 'A13_41';
        app.AlgorithmForCase135.Value = 'A135_41';
%         app.EnforcesymmetryaboutthexyplaneCheckBox.Value = true;
%         app.EnforcesymmetryaboutthexzplaneCheckBox.Value = true;
%         app.EnforcesymmetryabouttheyzplaneCheckBox.Value = true;
    case 'Barrel'
        app.EpsEditField.Value = 1e-10;
        app.ThicknessEditField.Value = 0.4;
        app.AlgorithmForCase13.Value = 'A13_31';
        app.AlgorithmForCase135.Value = 'A135_31';
        app.nurbsObjects{1,1}.nurbs = read_g2(['NURBSgeometries/g2files/' model '.g2']);
    case 'randomCubicHole'
        app.EpsEditField.Value = 1e-10;
        app.ThicknessEditField.Value = 1;
        app.nurbsObjects{1,1}.nurbs = read_g2(['NURBSgeometries/g2files/' model '.g2']);
        app.UseaveragenormalvectorsCheckBox.Value = false;
    case 'Mutter'
%         app.nurbsObjects{1,1}.nurbs = read_g2('../IGA-geometries/Mutter/Mutter.g2');
        app.nurbsObjects{1,1}.nurbs = getMutterData();
        app.Eps = 1e-4;
        app.EpsEditField.Value = app.Eps;

        app.ThicknessEditField.Value = 6;
        app.AnglethresholdEditField.Value = 140; % Threshold for a "sharp" angle
        app.PrioritizeleastnormaldeviationCheckBox.Value = false;
        app.MaintaincollapsednessCheckBox.Value = false; % Threshold for a "sharp" angle
        app.AlgorithmForCase13.Value = 'A13_41';
        app.AlgorithmForCase135.Value = 'A135_41';
        app.CaseswithsingularitiesEditField.Value = 0;
        app.ColorfunctionEditField.Value = '@(x) log10(abs(norm2(x(:,[1,3]))-6))';
    case 'BeTSSiM1'
        app.EpsEditField.Value = 1e-10;
        app.ThicknessEditField.Value = 1.5;
        app.AlgorithmForCase13.Value = 'A13_31';
        app.AlgorithmForCase135.Value = 'A135_31';
        app.nurbsObjects{1,1}.nurbs = read_g2(['NURBSgeometries/g2files/' model '.g2']);
    case 'BeTSSiM2'
        app.EpsEditField.Value = 1e-10;
        app.ThicknessEditField.Value = 1;
        app.AnglethresholdEditField.Value = 120; % Threshold for a "sharp" angle
%         app.ColorfunctionEditField.Value = '@(x) log10(abs(norm2(x(:,1:3))-1))';
%         app.AnglethresholdEditField.Value = 140; % Threshold for a "sharp" angle
%         app.AlgorithmForCase13.Value = 'A13_41';
%         app.AlgorithmForCase135.Value = 'A135_41';
        options.t = 2;
        app.nurbsObjects{1,1}.nurbs = getBeTSSiM2Data(options);
    case 'BeTSSiM3'
        app.EpsEditField.Value = 1e-10;
        app.ThicknessEditField.Value = 1.5;
        app.AnglethresholdEditField.Value = 140; % Threshold for a "sharp" angle
%         app.ColorfunctionEditField.Value = '@(x) log10(abs(norm2(x(:,1:3))-1))';
%         app.AnglethresholdEditField.Value = 140; % Threshold for a "sharp" angle
%         app.AlgorithmForCase13.Value = 'A13_41';
%         app.AlgorithmForCase135.Value = 'A135_41';
        app.nurbsObjects{1,1}.nurbs = read_g2(['NURBSgeometries/g2files/' model '.g2']);
        app.EnforcesymmetryaboutthexyplaneCheckBox.Value = true;
        app.EnforcesymmetryaboutthexzplaneCheckBox.Value = true;
        app.AlgorithmForCase13.Value = 'A13_41';
        app.UseaveragenormalvectorsCheckBox = false;
        view(app.UIAxes,[162,10]);
    case {'BeTSSiM4','BeTSSiM4_parm2'}
        if model(end) == '2'
            options.parm = 2;
            app.AnglethresholdEditField.Value = 150; % Threshold for a "sharp" angle
        else
            options.parm = 1;
            app.AlgorithmForCase135.Value = 'A135_31';
            app.AlgorithmForCase13.Value = 'A13_31';
            app.UseaveragenormalvectorsCheckBox.Value = false;
            app.AnglethresholdEditField.Value = 135; % Threshold for a "sharp" angle
        end
%         options.theta = 2*pi;r
%         options.theta_eta = 0.9*2*pi;
%         options.R = 2;
        options.t = 0.5;
        app.nurbsObjects{1,1}.nurbs = getBeTSSiM4Data(options);
%         app.nurbsObjects{1,1}.nurbs = read_g2(['NURBSgeometries/g2files/' model '.g2']);

        view(app.UIAxes,[135,25]);
        app.ThicknessEditField.Value = 1;
        app.AngledeviationEditField.Value = 40; % Threshold for a "sharp" angle
        app.ColorfunctionEditField.Value = '@(x) log10(abs(norm2(x(:,1:3)-[1,1,1]*0.5/2)-3))';
    case 'BCA'
        load('BeTSSi_BCA_p2_unRefined.mat','nurbs')
        app.nurbsObjects{1,1}.nurbs = nurbs;
        app.ThicknessEditField.Value = 2;
        app.MaintaincollapsednessCheckBox.Value = false; % Threshold for a "sharp" angle
        % sharpAngle = 100*pi/180; % Threshold for a "sharp" angle
        app.AnglethresholdEditField.Value = 120; % Threshold for a "sharp" angle
        % sharpAngle = 125*pi/180; % Threshold for a "sharp" angle
        % sharpAngle = 135*pi/180; % Threshold for a "sharp" angle
        % sharpAngle = 160*pi/180; % Threshold for a "sharp" angle
        % sharpAngle = 161*pi/180; % Threshold for a "sharp" angle % In order to include connections over depth rudders
        % sharpAngle = 173*pi/180; % Threshold for a "sharp" angle
        app.EnforcesymmetryaboutthexzplaneCheckBox.Value = true;
        app.UseaveragenormalvectorsCheckBox.Value = false;
        view(app.UIAxes,[162,10]);
    case 'Bridge'
        nurbs = read_g2('../../OneDrive - SINTEF/SINTEF/DT/bridge-quadratic/bridge-quadratic.g2');
        nurbs = translateNURBS(nurbs,[-55,0,0]);
        app.nurbsObjects{1,1}.nurbs = nurbs;
        app.EnforcesymmetryabouttheyzplaneCheckBox.Value = true;
%         app.nurbsObjects{1,1}.nurbs = read_g2('../../OneDrive/SINTEF/DT/bridge-quadratic/bridge-quadratic.g2');
        app.AlgorithmForCase135.Value = 'A135_11';
        app.AlgorithmForCase13.Value = 'A13_11';
        app.xminEditField.Value = -150;
        app.xmaxEditField.Value = 150;
        app.yminEditField.Value = -100;
        app.ymaxEditField.Value = 200;
        app.zminEditField.Value = 0;
        app.zmaxEditField.Value = 100;
        app.Seth_maxEditField.Value = 1;
%         app.nurbsObjects{1,1}.nurbs = app.nurbsObjects{1,1}.nurbs([93,94,120,121,127,128,141,142,146,147]);
%         app.nurbsObjects{1,1}.nurbs = app.nurbsObjects{1,1}.nurbs([1,2,28,29,35,36,50,55,49,54,114]);
%         app.nurbsObjects{1,1}.nurbs = app.nurbsObjects{1,1}.nurbs([1:3,26:65]);
%         app.Eps = 1e-4;
%         app.EpsEditField.Value = app.Eps;
        app.UseaveragenormalvectorsCheckBox.Value = false;
        app.MaintaincollapsednessCheckBox.Value = false; % Threshold for a "sharp" angle

        app.ThicknessEditField.Value = 54.25;
        app.AnglethresholdEditField.Value = 120; % Threshold for a "sharp" angle
    case 'Capablanca'
        load('../../results/ASIGA/Capablanca/nurbs.mat','nurbs')
        app.nurbsObjects{1,1}.nurbs = nurbs;

        app.ThicknessEditField.Value = 1.1;
        app.AnglethresholdEditField.Value = 120; % Threshold for a "sharp" angle
    case 'Capablanca_woRooks'
        load('../../results/ASIGA/Capablanca/nurbs_woRook.mat','nurbs')
        app.nurbsObjects{1,1}.nurbs = nurbs;
        app.ThicknessEditField.Value = 1.1;
        app.AnglethresholdEditField.Value = 120; % Threshold for a "sharp" angle
    case 'Pawn'
        app.nurbsObjects{1,1}.nurbs = read_g2('../IGA-geometries/Chess/Pawn.g2');
        app.Eps = 1e-4;
        app.EpsEditField.Value = app.Eps;

        app.ThicknessEditField.Value = 1.1;
        app.AnglethresholdEditField.Value = 120; % Threshold for a "sharp" angle
    case 'S1_interior'
        app.nurbsObjects{1,1}.nurbs = getEllipsoidData('parm',2);
        app.AlgorithmForCase1.Value = 'A1_31';
        app.AlgorithmForCase13.Value = 'A13_41';
        app.AlgorithmForCase135.Value = 'A135_41';
        app.AlgorithmForCase135.Value = 'A135_41';
        app.nurbsObjects{1,1}.nurbs = flipNURBSparametrization(app.nurbsObjects{1,1}.nurbs,1);
        app.ThicknessEditField.Value = 0.5;
        app.AlphavalueSlider.Value = 0.5;
        app.AnglethresholdEditField.Value = 140; % Threshold for a "sharp" angle
        app.ColorfunctionEditField.Value = '@(x) log10(abs(norm2(x(:,1:3))-1))';
        
    otherwise
        app.EpsEditField.Value = 1e-10;
        if numel(model) >= 6 && strcmp(model(1:6), 'BeTSSi')
            app.ThicknessEditField.Value = 1.5;
        else
            app.ThicknessEditField.Value = 0.1;
        end
        app.ColorfunctionEditField.Value = '@(x) log10(abs(norm2(x(:,1:3))-1))';
        app.nurbsObjects{1,1}.nurbs = read_g2(['NURBSgeometries/g2files/' model '.g2']);
        app.AnglethresholdEditField.Value = 140; % Threshold for a "sharp" angle
end