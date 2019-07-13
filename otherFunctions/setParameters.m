function parms = setParameters(model,parm)

switch model
    case 'SS_P'
        setPulsatingSphereParameters
    case 'EL'
        setEllipsoidParameters   
    case 'SS'
        setSSParameters
    case {'S1','S1_P','S1_P2'}
        setS1Parameters
    case 'S3'
        setS3Parameters
    case 'Torus'
        setTorusParameters
    case {'Cube','Cube_P'}
        setCubeParameters
    case 'S5'
        setS5Parameters
    case 'IL'
        setIhlenburgParameters
    case {'M1','M1_P'}
        setBeTSSi_M1Parameters
    case {'M2','M2_P'}
        setBeTSSi_M2Parameters
    case {'M3','M3_P'}
        setBeTSSi_M3Parameters
    case {'PH','PH_P'}
        setBeTSSi_PHParameters
    case {'M4','M4_P'}
        setBeTSSi_M4Parameters
    case {'M5A','M5B','M5A_P','M5B_P'}
        setBeTSSi_M5Parameters
    case {'MS','MS_P'}
        setMockShellParameters
        if ~isempty(parm)
            mult = parm;
        end
        L = mult*R_o*pi/2; 
    case 'TAP'
        setTAPParameters
    case {'Barrel','Barrel_P'}
        setBarrelParameters
    case {'BC','BC_P','BCA','BCA_P'}
        setBCParameters
end
variables = whos;
for i = 1:length(variables)
    parms.(variables(i).name) = eval(variables(i).name);
end