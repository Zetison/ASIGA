% Create g2 files from available matlab geometries
close all
model = {'S1','S1patched','M1','M2','M3','M4','M5A','M5B','PH','Torus','Cube','Shirron','Barrel','BCA'};
model = {'Cube'};

for i = 1:numel(model)
    switch model{i}
        case 'S1'
            varCol = setS1Parameters;
            t = varCol{1}.R_i - varCol{2}.R_i;
            nurbs = getEllipsoidData();
        case 'S1patched'
            varCol = setS1Parameters;
            t = varCol{1}.R_i - varCol{2}.R_i;
            nurbs = getEllipsoidData('parm',2); 
        case 'M1'
            nurbs = getBeTSSiM1Data();
        case 'M2'
            nurbs = getBeTSSiM2Data();
        case 'M3'
            nurbs = getBeTSSiM3Data();
        case 'M4'
            nurbs = getBeTSSiM4Data();
        case 'M5A'
            nurbs = getBeTSSiM5Data('type','A');
        case 'M5B'
            nurbs = getBeTSSiM5Data('type','B');
        case 'PH'
            nurbs = getBeTSSiPHData();
        case 'Torus'
            nurbs = getTorusData();
        case 'Cube'
            nurbs = getCubeData(); 
        case 'Shirron'
            varCol = setShirronParameters;
            nurbs = getBeTSSiM3Data('R1',varCol{1}.R_o,'R2',varCol{1}.R_o, 'L', varCol{1}.L);
        case 'Barrel'
            nurbs = getBarrelData();
        case 'BCA'
            for degree = 2:15
                load(['NURBSgeometries/BCAdata/BeTSSi_BCA_p' num2str(degree)])
                write_g2(nurbs,['NURBSgeometries/g2files/' model{i} '_p' num2str(degree)])
            end
    end
    if ~strcmp(model{i},'BCA')
        write_g2(nurbs,['NURBSgeometries/g2files/' model{i}])
    end
end
plotNURBS(nurbs)
camlight