function app = getS2Vmodel(app,model)
app.AlgorithmForCase13.Value = 'A13_21';
app.AlgorithmForCase135.Value = 'A135_21';
switch model
    case 'Mutter'
        app.nurbsObjects{1,1}.nurbs = read_g2('../IGA-geometries/Mutter/Mutter.g2');
        app.Eps = 1e-4;
        app.EpsEditField.Value = app.Eps;

        app.ThicknessEditField.Value = 1.1;
        app.AnglethresholdEditField.Value = 120; % Threshold for a "sharp" angle
        app.ColorfunctionEditField.Value = '@(x) log10(abs(norm2(x(:,[1,3]))-6))';
    case 'Bridge'
        app.nurbsObjects{1,1}.nurbs = read_g2('../../OneDrive/SINTEF/DT/bridge-quadratic/bridge-quadratic.g2');
%         app.Eps = 1e-4;
%         app.EpsEditField.Value = app.Eps;

        app.ThicknessEditField.Value = 1.1;
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
    case 'BCA'
        load('BeTSSi_BCA_p2_unRefined.mat','nurbs')
        app.nurbsObjects{1,1}.nurbs = nurbs;
        app.ThicknessEditField.Value = 1.1;
        % sharpAngle = 100*pi/180; % Threshold for a "sharp" angle
        app.AnglethresholdEditField.Value = 120; % Threshold for a "sharp" angle
        % sharpAngle = 125*pi/180; % Threshold for a "sharp" angle
        % sharpAngle = 135*pi/180; % Threshold for a "sharp" angle
        % sharpAngle = 160*pi/180; % Threshold for a "sharp" angle
        % sharpAngle = 161*pi/180; % Threshold for a "sharp" angle % In order to include connections over depth rudders
        % sharpAngle = 173*pi/180; % Threshold for a "sharp" angle
    case 'S1_interior'
        app.nurbsObjects{1,1}.nurbs = getEllipsoidData('parm',2);
        app.nurbsObjects{1,1}.nurbs = flipNURBSparametrization(app.nurbsObjects{1,1}.nurbs,1);
        app.AnglethresholdEditField.Value = 140; % Threshold for a "sharp" angle
        app.ColorfunctionEditField.Value = '@(x) log10(abs(norm2(x(:,1:3))-1))';
    case 'S1'
        app.EpsEditField.Value = 1e-10;
        app.ThicknessEditField.Value = 0.1;
        app.ColorfunctionEditField.Value = '@(x) log10(abs(norm2(x(:,1:3))-1))';
        app.nurbsObjects{1,1}.nurbs = read_g2(['NURBSgeometries/g2files/' model '.g2']);
        app.AnglethresholdEditField.Value = 140; % Threshold for a "sharp" angle
        app.AlgorithmForCase13.Value = 'A13_41';
        app.AlgorithmForCase135.Value = 'A135_41';
    case 'BeTSSiM1'
        app.EpsEditField.Value = 1e-10;
        app.ThicknessEditField.Value = 1.5;
        app.AlgorithmForCase13.Value = 'A13_31';
        app.AlgorithmForCase135.Value = 'A135_31';
        app.nurbsObjects{1,1}.nurbs = read_g2(['NURBSgeometries/g2files/' model '.g2']);
    case 'BeTSSiM2'
        app.EpsEditField.Value = 1e-10;
        app.ThicknessEditField.Value = 0.1;
%         app.ColorfunctionEditField.Value = '@(x) log10(abs(norm2(x(:,1:3))-1))';
%         app.AnglethresholdEditField.Value = 140; % Threshold for a "sharp" angle
%         app.AlgorithmForCase13.Value = 'A13_41';
%         app.AlgorithmForCase135.Value = 'A135_41';
        options.t = 2;
        app.nurbsObjects{1,1}.nurbs = getBeTSSiM2Data(options);
    case {'BeTSSiM4','BeTSSiM4_parm2'}
        if model(end) == '2'
            options.parm = 2;
        else
            options.parm = 1;
            app.AlgorithmForCase135.Value = 'A135_31';
            app.AlgorithmForCase13.Value = 'A13_31';
        end
%         options.theta = 2*pi;r
%         options.theta_eta = 0.9*2*pi;
%         options.R = 2;
        options.t = 0.5;
        app.nurbsObjects{1,1}.nurbs = getBeTSSiM4Data(options);
%         app.nurbsObjects{1,1}.nurbs = read_g2(['NURBSgeometries/g2files/' model '.g2']);

        view(app.UIAxes,[135,25]);
        app.ThicknessEditField.Value = 1;
        app.AnglethresholdEditField.Value = 150; % Threshold for a "sharp" angle
        app.AngledeviationEditField.Value = 45; % Threshold for a "sharp" angle
        app.ColorfunctionEditField.Value = '@(x) log10(abs(norm2(x(:,1:3)-[1,1,1]*0.5/2)-3))';
    case 'randomCubicHole'
        app.EpsEditField.Value = 1e-10;
        app.ThicknessEditField.Value = 1;
        app.nurbsObjects{1,1}.nurbs = read_g2(['NURBSgeometries/g2files/' model '.g2']);
    case 'Barrel'
        app.EpsEditField.Value = 1e-10;
        app.ThicknessEditField.Value = 0.4;
        app.AlgorithmForCase135.Value = 'A135_31';
        app.nurbsObjects{1,1}.nurbs = read_g2(['NURBSgeometries/g2files/' model '.g2']);
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