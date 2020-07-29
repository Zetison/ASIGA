function viewAsEl = getView(type)
if nargin < 1
    type = 0;
end
switch type
    case 0
        viewAsEl = [18,10];
    case 1
        viewAsEl = [0,90]; % from above
    case 2
        viewAsEl = [11,25]; % BeTSSi M2
    case 3
        viewAsEl = [46,32]; % BeTSSi M4
    case 4
        viewAsEl = [106,26]; % sphere and cube
    case 5
        viewAsEl = [18+180,10]; % BeTSSi M1
end