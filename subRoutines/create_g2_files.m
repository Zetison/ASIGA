% Create g2 files from available matlab geometries

model = {'S1','S1patched','M1','M2','M3','M4','M5A','M5B','PH','Torus','Cube','Shirron','Barrel','BCA'};
model = {'BCA'};

for i = 1:numel(model)
    switch model{i}
        case 'S1'
            setS1Parameters
            nurbs = getEllipsoidalShellData(R_o,R_o,R_o,t,'Zaxis');
        case 'S1patched'
            setS1Parameters
            nurbs = getSphericalShellDataPatched(R_o, t); 
        case 'M1'
            setBeTSSi_M1Parameters
            nurbs = getModel1Data(R_i, R_o, L, 1/3, 2/3);
        case 'M2'
            setBeTSSi_M2Parameters
            nurbs = getModel2Data(R_o,t,L,theta2);
        case 'M3'
            setBeTSSi_M3Parameters
            nurbs = getModel3Data(R_o1, R_o2, t, L);
        case 'M4'
            setBeTSSi_M4Parameters
            nurbs = getModel4Data(R_o, t, 1/3, 2/3);
        case 'M5A'
            setBeTSSi_M5Parameters
            nurbs = getModel5Data(R_o, 1/3, 2/3, L, l, 'Xaxis');
        case 'M5B'
            setBeTSSi_M5Parameters
            nurbs = getModel5Data(R_o, 1/3, 2/3, L, l, 'Zaxis');
        case 'PH'
            setBeTSSi_PHParameters
            nurbs = getPHData(R_1, R_2, t, gd);
        case 'Torus'
            setTorusParameters
            nurbs = getTorusData(r_o,r_i);
        case 'Cube'
            setCubeParameters
            nurbs = getCubeData(a); 
        case 'Shirron'
            setShirronParameters
            nurbs = getModel3Data(R_o, R_o, t, L, 4);
        case 'Barrel'
            setBarrelParameters
            nurbs = getBarrelData(R_o,R_i,1/3,2/3,L);
        case 'BCA'
            setBCParameters
            for degree = 2:15
                load(['NURBSgeometries/BCAdata/BeTSSi_BCA_p' num2str(degree)])
                write_g2(nurbs,['NURBSgeometries/g2files/' model{i} '_p' num2str(degree)])
            end
    end
    if ~strcmp(model{i},'BCA')
        write_g2(nurbs,['NURBSgeometries/g2files/' model{i}])
    end
end