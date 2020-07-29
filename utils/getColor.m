function color = getColor(type)

switch type
    case 0
        color = getColor(1);
    case 1
        color = [208 213 219]/255; % polished aluminum
    case 2
        color = [152 215 112]/255; % Cottrell2006iat green 1
    case 3
        color = [162 211 78]/255;  % Cottrell2006iat green 2
    case 4
        color = [122 127 128]/255; % metal gray
    case 5
        color = [132 135 137]/255; % aluminum
    case 6
        color = [0,80,158]/255; % NTNU logo color
    case 7
        color = [0, 60, 101]/255; % SINTEF logo color
    case 8
        color = [0, 68, 123]/255; % FFI logo color
    case 9
        color = [0, 77, 145]/255; % background FFI color
end