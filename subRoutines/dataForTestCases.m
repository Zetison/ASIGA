counter = 1;
testCases(counter).solid_degreeElevArr = [0 0 1];
testCases(counter).fluid_degreeElevArr = [0 0 0];
newXiKnots = 0; % 3 -> 14
newEtaKnots = 0; % 30 -> 14
testCases(counter).fluid_newKnotsArr = [newXiKnots newEtaKnots 2]; % - - 2
testCases(counter).solid_newKnotsArr = [newXiKnots newEtaKnots 2];
testCases(counter).noRadShapeFuncInfElements = 3;
testCases(counter).R_a = 6;
testCases(counter).infiniteElementFormulation = 'PGC'; 
counter = counter + 1;
return

if analyzeFEMcase
    counter = 1;
    %% polynomial enrichment in fluid
    testCases(counter).solid_degreeElevArr = [0 0 1];
    testCases(counter).fluid_degreeElevArr = [0 0 0];
    newXiKnots = 2^(M-1)-1; % 3
    newEtaKnots = 2^(M-1)-1; % 30
    testCases(counter).fluid_newKnotsArr = [newXiKnots newEtaKnots 0]; % - - 2
    testCases(counter).solid_newKnotsArr = [newXiKnots newEtaKnots 0];
    testCases(counter).noRadShapeFuncInfElements = 3;
    testCases(counter).R_a = 6;
    testCases(counter).infiniteElementFormulation = 'PGU'; 
    counter = counter + 1;

    testCases(counter) = testCases(counter-1);
    testCases(counter).fluid_degreeElevArr = [0 0 1];
    counter = counter + 1;

    testCases(counter) = testCases(counter-2);
    testCases(counter).fluid_degreeElevArr = [0 0 2];
    counter = counter + 1;

    %% polynomial enrichment in solid
    testCases(counter).solid_degreeElevArr = [0 0 0];
    testCases(counter).fluid_degreeElevArr = [0 0 1];
    newXiKnots = 2^(M-1)-1; % 3
    newEtaKnots = 2^(M-1)-1; % 30
    testCases(counter).fluid_newKnotsArr = [newXiKnots newEtaKnots 0]; % - - 2
    testCases(counter).solid_newKnotsArr = [newXiKnots newEtaKnots 0];
    testCases(counter).noRadShapeFuncInfElements = 3;
    testCases(counter).R_a = 6;
    testCases(counter).infiniteElementFormulation = 'PGU'; 
    counter = counter + 1;


    testCases(counter) = testCases(counter-1);
    testCases(counter).solid_degreeElevArr = [0 0 1];
    counter = counter + 1;

    testCases(counter) = testCases(counter-2);
    testCases(counter).solid_degreeElevArr = [0 0 2];
    counter = counter + 1;


    %% Infinite element basis function enrichment
    testCases(counter).solid_degreeElevArr = [0 0 1];
    testCases(counter).fluid_degreeElevArr = [0 0 1];
    newXiKnots = 2^(M-1)-1; % 3
    newEtaKnots = 2^(M-1)-1; % 30
    testCases(counter).fluid_newKnotsArr = [newXiKnots newEtaKnots 0]; % - - 2
    testCases(counter).solid_newKnotsArr = [newXiKnots newEtaKnots 0];
    testCases(counter).noRadShapeFuncInfElements = 1;
    testCases(counter).R_a = 6;
    testCases(counter).infiniteElementFormulation = 'PGU'; 
    counter = counter + 1;


    testCases(counter) = testCases(counter-1);
    testCases(counter).noRadShapeFuncInfElements = 3;
    counter = counter + 1;

    testCases(counter) = testCases(counter-2);
    testCases(counter).noRadShapeFuncInfElements = 6;
    counter = counter + 1;


    %% Infinite element basis function enrichment
    testCases(counter).solid_degreeElevArr = [0 0 1];
    testCases(counter).fluid_degreeElevArr = [0 0 1];
    newXiKnots = 2^(M-1)-1; % 3
    newEtaKnots = 2^(M-1)-1; % 30
    testCases(counter).fluid_newKnotsArr = [newXiKnots newEtaKnots 0]; % - - 2
    testCases(counter).solid_newKnotsArr = [newXiKnots newEtaKnots 0];
    testCases(counter).noRadShapeFuncInfElements = 3;
    testCases(counter).R_a = 6;
    testCases(counter).infiniteElementFormulation = 'PGC'; 
    counter = counter + 1;

    % infiniteElementFormulations = {'PGC', 'PGU', 'BGU', 'BGC'};

    testCases(counter) = testCases(counter-1);
    testCases(counter).infiniteElementFormulation = 'PGU';
    counter = counter + 1;

    testCases(counter) = testCases(counter-2);
    testCases(counter).infiniteElementFormulation = 'BGU';
    counter = counter + 1;

    testCases(counter) = testCases(counter-3);
    testCases(counter).infiniteElementFormulation = 'BGC';
    counter = counter + 1;
else
    counter = 1;
    %% polynomial enrichment in fluid
    testCases(counter).solid_degreeElevArr = [0 0 1];
    testCases(counter).fluid_degreeElevArr = [0 0 0];
    newXiKnots = 2*2^(M-1)-1; % 3 -> 14
    newEtaKnots = 2*2^(M-1)-1; % 30 -> 14
    testCases(counter).fluid_newKnotsArr = [newXiKnots newEtaKnots 0]; % - - 2
    testCases(counter).solid_newKnotsArr = [newXiKnots newEtaKnots 0];
    testCases(counter).noRadShapeFuncInfElements = 3;
    testCases(counter).R_a = 6;
    testCases(counter).infiniteElementFormulation = 'PGU'; 
    counter = counter + 1;

    testCases(counter) = testCases(counter-1);
    testCases(counter).fluid_degreeElevArr = [0 0 1];
    counter = counter + 1;

    testCases(counter) = testCases(counter-2);
    testCases(counter).fluid_degreeElevArr = [0 0 2];
    counter = counter + 1;

    %% polynomial enrichment in solid
    testCases(counter).solid_degreeElevArr = [0 0 0];
    testCases(counter).fluid_degreeElevArr = [0 0 1];
    newXiKnots = 2*2^(M-1)-1; % 3
    newEtaKnots = 2*2^(M-1)-1; % 30
    testCases(counter).fluid_newKnotsArr = [newXiKnots newEtaKnots 0]; % - - 2
    testCases(counter).solid_newKnotsArr = [newXiKnots newEtaKnots 0];
    testCases(counter).noRadShapeFuncInfElements = 3;
    testCases(counter).R_a = 6;
    testCases(counter).infiniteElementFormulation = 'PGU'; 
    counter = counter + 1;


    testCases(counter) = testCases(counter-1);
    testCases(counter).solid_degreeElevArr = [0 0 1];
    counter = counter + 1;

    testCases(counter) = testCases(counter-2);
    testCases(counter).solid_degreeElevArr = [0 0 2];
    counter = counter + 1;


    %% Infinite element basis function enrichment
    testCases(counter).solid_degreeElevArr = [0 0 1];
    testCases(counter).fluid_degreeElevArr = [0 0 1];
    newXiKnots = 2*2^(M-1)-1; % 3
    newEtaKnots = 2*2^(M-1)-1; % 30
    testCases(counter).fluid_newKnotsArr = [newXiKnots newEtaKnots 0]; % - - 2
    testCases(counter).solid_newKnotsArr = [newXiKnots newEtaKnots 0];
    testCases(counter).noRadShapeFuncInfElements = 1;
    testCases(counter).R_a = 6;
    testCases(counter).infiniteElementFormulation = 'PGU'; 
    counter = counter + 1;


    testCases(counter) = testCases(counter-1);
    testCases(counter).noRadShapeFuncInfElements = 3;
    counter = counter + 1;

    testCases(counter) = testCases(counter-2);
    testCases(counter).noRadShapeFuncInfElements = 6;
    counter = counter + 1;


    %% Infinite element basis function enrichment
    testCases(counter).solid_degreeElevArr = [0 0 1];
    testCases(counter).fluid_degreeElevArr = [0 0 1];
    newXiKnots = 2*2^(M-1)-1; % 3
    newEtaKnots = 2*2^(M-1)-1; % 30
    testCases(counter).fluid_newKnotsArr = [newXiKnots newEtaKnots 0]; % - - 2
    testCases(counter).solid_newKnotsArr = [newXiKnots newEtaKnots 0];
    testCases(counter).noRadShapeFuncInfElements = 3;
    testCases(counter).R_a = 6;
    testCases(counter).infiniteElementFormulation = 'PGC'; 
    counter = counter + 1;

    % infiniteElementFormulations = {'PGC', 'PGU', 'BGU', 'BGC'};

    testCases(counter) = testCases(counter-1);
    testCases(counter).infiniteElementFormulation = 'PGU';
    counter = counter + 1;

    testCases(counter) = testCases(counter-2);
    testCases(counter).infiniteElementFormulation = 'BGU';
    counter = counter + 1;

    testCases(counter) = testCases(counter-3);
    testCases(counter).infiniteElementFormulation = 'BGC';
    counter = counter + 1;
end


