function app = getS2Vmodel(app,model)
switch model
    case 'Mutter'
        app.nurbsObjects{1,1}.nurbs = read_g2('../IGA-geometries/Mutter/Mutter.g2');
        app.Eps = 1e-4;
        app.EpsEditField.Value = app.Eps;

        app.ThicknessEditField.Value = 1.1;
        app.AnglethresholdEditField.Value = 120; % Threshold for a "sharp" angle
        app.ColorfunctionEditField.Value = '@(x) log10(abs(norm2(x(:,[1,3]))-6))';
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
        app.nurbsObjects{1,1}.nurbs = getEllipsoidData('parm',2,'S2V_algorithm',{'A1_2'});
        app.nurbsObjects{1,1}.nurbs = flipNURBSparametrization(app.nurbsObjects{1,1}.nurbs,1);
        app.AnglethresholdEditField.Value = 140; % Threshold for a "sharp" angle
        app.ColorfunctionEditField.Value = '@(x) log10(abs(norm2(x(:,1:3))-1))';
    case 'BeTSSiM4_parm2'
        options.parm = 2;
%         options.theta = 2*pi;r
%         options.theta_eta = 0.9*2*pi;
%         options.R = 2;
        options.t = 0.5;
        app.nurbsObjects{1,1}.nurbs = getBeTSSiM4Data(options);
%         app.nurbsObjects{1,1}.nurbs = read_g2(['NURBSgeometries/g2files/' model '.g2']);
        app.ThicknessEditField.Value = 0.1;
        app.AnglethresholdEditField.Value = 140; % Threshold for a "sharp" angle
        app.AngledeviationEditField.Value = 45; % Threshold for a "sharp" angle
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