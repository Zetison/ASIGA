

scatteringCase = 'MS'; % 'BI' = Bistatic scattering, 'MS' = Monostatic scattering

model = 'TAP'; % The TAP model (Mock shell) described in Richard Hodges book "Underwater acoustics: analysis, design and performance of sonar"

calculateFarFieldPattern = 1;

f = 299;             % Frequency
f = 100;             % Frequency

alpha = (0:0.1:180)*pi/180;
% alpha = (0:180)*pi/180;
beta = 0;   

plot2Dgeometry = 0;
plot3Dgeometry = 0;

%% Plot result in paraview
if 0
    if true
        coreMethod = 'IGA';
        M = 4;
        degree = 2:3;
        loopParameters = {'M','degree','BC'};
    else
        coreMethod = 'SEM';
        M = 1;
        degree = 2:5;
        loopParameters = {'degree','BC'};
    end
    N = 4;
    formulation = 'BGU';
    method = 'IE';
else
    M = 1:4;
    parm = linspace(0.1,3,10); % max ~5
    parm = 0.2;
    method = 'MFS';
    formulation = '';
    degree = 2;
    loopParameters = {'M','BC','parm'};
    computeCondNumber = 0;
end
BC = {'SHBC'};

collectIntoTasks
