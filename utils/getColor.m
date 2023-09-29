function color = getColor(type)
if nargin < 1
    type = 1;
end
switch type
    case 0
        color = getColor(1);
    case 1
        color = [208 213 219]; % polished aluminum
    case 2
        color = [152 215 112]; % Cottrell2006iat green 1
    case 3
        color = [162 211 78];  % Cottrell2006iat green 2
    case 4
        color = [122 127 128]; % metal gray
    case 5
        color = [132 135 137]; % aluminum
    case 6
        color = [0,80,158]; % NTNU logo color
    case 7
        color = [0, 60, 101]; % SINTEF logo color
    case 8
        color = [0, 68, 123]; % FFI logo color
    case 9
        color = [0, 77, 145]; % background FFI color
    case 10
        color = [173, 216, 230]; % water
    case 11 
        color = [0,127,0]; % PML color
    case 12 
        color = [255,0,0]*0.3; % myRed
    case 13
        color = [0,255,0]*0.5; % myGreen
end
color = color/255;