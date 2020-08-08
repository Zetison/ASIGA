function [adjacentElements,xi_x_tArr,eta_x_tArr] = getAdjacentElements(e_x,xi_x,eta_x,Xi_e_x,Eta_e_x,eNeighbour,Eps)
% find adjecent elements around the colocation point x in element e_x
% (located at (xi_x,eta_x) in the parameter space)
adjacentElements = NaN(1,8);
adjacentElements(1) = e_x;
xi_x_tArr = NaN(1,8);
eta_x_tArr = NaN(1,8);
xi_x_tArr(1) = parametric2parentSpace(Xi_e_x, xi_x);
eta_x_tArr(1) = parametric2parentSpace(Eta_e_x, eta_x);

pointIsInEast = false;
pointIsInNorth = false;
pointIsInWest = false;
pointIsInSouth = false;
if abs(xi_x_tArr(1) - (-1)) < Eps
    pointIsInWest = true;
elseif abs(xi_x_tArr(1) - 1) < Eps
    pointIsInEast = true;
end
if abs(eta_x_tArr(1) - (-1)) < Eps
    pointIsInSouth = true;
elseif abs(eta_x_tArr(1) - 1) < Eps
    pointIsInNorth = true;
end
cornerPoints = [ 1, 1;
                -1, 1;
                -1,-1;
                 1,-1];
counter = 2;
if pointIsInEast && pointIsInNorth % top right corner (xi_x_t = 1, eta_x_t = 1)
    counter3 = 1;
    e_prev = e_x;
    e_next = eNeighbour(e_x,1,1); % element to the east
    if isnan(e_next)
        e_next = eNeighbour(e_x,4,1); % element to the south
        counter3 = counter3 - 1;
    end
    while e_next ~= e_x % go counter clockwise around vertex to end up back at e_x
        counter2 = 3; % index of element to the west
        e_temp = eNeighbour(e_next,counter2,1);
        while e_temp ~= e_prev
            counter2 = counter2 + 1;
            e_temp = eNeighbour(e_next,mod(counter2-1,4)+1,1);
        end

        adjacentElements(counter) = e_next;
        twist = mod(4-eNeighbour(e_next,mod(counter2-1,4)+1,2),4);
        xi_vec = cornerPoints(mod(counter3+1-twist-1,4)+1,:);
        xi_x_tArr(counter)  = xi_vec(1);
        eta_x_tArr(counter) =  xi_vec(2);
        
        e_prev = e_next;
        e_next = eNeighbour(e_prev,mod(counter2+3-1,4)+1,1);
        if isnan(e_next)
            e_next = eNeighbour(e_prev,mod(counter2+2-1,4)+1,1);
            counter3 = counter3 - 1;
        end
        counter = counter + 1;
        counter3 = counter3 + 1 - twist;
    end
elseif pointIsInWest && pointIsInNorth % top left corner (xi_x_t = -1, eta_x_t = 1)
    counter3 = 2;
    e_prev = e_x;
    e_next = eNeighbour(e_x,2,1); % element to the north
    if isnan(e_next)
        e_next = eNeighbour(e_x,1,1); % element to the east
        counter3 = counter3 - 1;
    end
    while e_next ~= e_x % go counter clockwise around vertex to end up back at e_x
        counter2 = 4; % index of element to the south
        e_temp = eNeighbour(e_next,counter2,1);
        while e_temp ~= e_prev
            counter2 = counter2 + 1;
            e_temp = eNeighbour(e_next,mod(counter2-1,4)+1,1);
        end

        adjacentElements(counter) = e_next;
        twist = mod(4-eNeighbour(e_next,mod(counter2-1,4)+1,2),4);
        xi_vec = cornerPoints(mod(counter3+1-twist-1,4)+1,:);
        xi_x_tArr(counter)  = xi_vec(1);
        eta_x_tArr(counter) =  xi_vec(2);
        
        e_prev = e_next;
        e_next = eNeighbour(e_prev,mod(counter2+3-1,4)+1,1);
        if isnan(e_next)
            e_next = eNeighbour(e_prev,mod(counter2+2-1,4)+1,1);
            counter3 = counter3 - 1;
        end
        counter = counter + 1;
        counter3 = counter3 + 1 - twist;
    end     
elseif pointIsInWest && pointIsInSouth % bottom left corner (xi_x_t = -1, eta_x_t = -1)
    counter3 = 3;
    e_prev = e_x;
    e_next = eNeighbour(e_x,3,1); % element to the west
    if isnan(e_next)
        e_next = eNeighbour(e_x,2,1); % element to the north
        counter3 = counter3 - 1;
    end
    while e_next ~= e_x % go counter clockwise around vertex to end up back at e_x
        counter2 = 1; % index of element to the east
        e_temp = eNeighbour(e_next,counter2,1);
        while e_temp ~= e_prev
            counter2 = counter2 + 1;
            e_temp = eNeighbour(e_next,mod(counter2-1,4)+1,1);
        end

        adjacentElements(counter) = e_next;
        twist = mod(4-eNeighbour(e_next,mod(counter2-1,4)+1,2),4);
        xi_vec = cornerPoints(mod(counter3+1-twist-1,4)+1,:);
        xi_x_tArr(counter)  = xi_vec(1);
        eta_x_tArr(counter) =  xi_vec(2);
        
        e_prev = e_next;
        e_next = eNeighbour(e_prev,mod(counter2+3-1,4)+1,1);
        if isnan(e_next)
            e_next = eNeighbour(e_prev,mod(counter2+2-1,4)+1,1);
            counter3 = counter3 - 1;
        end
        counter = counter + 1;
        counter3 = counter3 + 1 - twist;
    end     
elseif pointIsInEast && pointIsInSouth % bottom right corner (xi_x_t = 1, eta_x_t = -1)
    counter3 = 4;
    e_prev = e_x;
    e_next = eNeighbour(e_x,4,1); % element to the south
    if isnan(e_next)
        e_next = eNeighbour(e_x,3,1); % element to the west
        counter3 = counter3 - 1;
    end
    while e_next ~= e_x % go counter clockwise around vertex to end up back at e_x
        counter2 = 2; % index of element to the north
        e_temp = eNeighbour(e_next,counter2,1);
        while e_temp ~= e_prev
            counter2 = counter2 + 1;
            e_temp = eNeighbour(e_next,mod(counter2-1,4)+1,1);
        end

        adjacentElements(counter) = e_next;
        twist = mod(4-eNeighbour(e_next,mod(counter2-1,4)+1,2),4);
        xi_vec = cornerPoints(mod(counter3+1-twist-1,4)+1,:);
        xi_x_tArr(counter)  = xi_vec(1);
        eta_x_tArr(counter) =  xi_vec(2);
        
        e_prev = e_next;
        e_next = eNeighbour(e_prev,mod(counter2+3-1,4)+1,1);
        if isnan(e_next)
            e_next = eNeighbour(e_prev,mod(counter2+2-1,4)+1,1);
            counter3 = counter3 - 1;
        end
        counter = counter + 1;
        counter3 = counter3 + 1 - twist;
    end     
elseif pointIsInEast % xi_x_t = 1
    adjacentElements(2) = eNeighbour(e_x,1,1);
    twist = eNeighbour(e_x,1,2);
    switch twist
        case 0
            xi_x_tArr(2)  = -1;
            eta_x_tArr(2) =  eta_x_tArr(1);
        case 1
            xi_x_tArr(2)  =  eta_x_tArr(1);
            eta_x_tArr(2) =  1;
        case 2
            xi_x_tArr(2)  =  1;
            eta_x_tArr(2) = -eta_x_tArr(1);
        case 3
            xi_x_tArr(2)  = -eta_x_tArr(1);
            eta_x_tArr(2) = -1;
    end
elseif pointIsInNorth % eta_x_t = 1
    adjacentElements(2) = eNeighbour(e_x,2,1);
    twist = eNeighbour(e_x,2,2);
    switch twist
        case 0
            xi_x_tArr(2)  =  xi_x_tArr(1);
            eta_x_tArr(2) = -1;
        case 1
            xi_x_tArr(2)  = -1;
            eta_x_tArr(2) = -xi_x_tArr(1);
        case 2
            xi_x_tArr(2)  = -xi_x_tArr(1);
            eta_x_tArr(2) =  1;
        case 3
            xi_x_tArr(2)  =  1;
            eta_x_tArr(2) =  xi_x_tArr(1);
    end
elseif pointIsInWest % xi_x_t = -1
    adjacentElements(2) = eNeighbour(e_x,3,1);
    twist = eNeighbour(e_x,3,2);
    switch twist
        case 0
            xi_x_tArr(2)  =  1;
            eta_x_tArr(2) =  eta_x_tArr(1);
        case 1
            xi_x_tArr(2)  =  eta_x_tArr(1);
            eta_x_tArr(2) = -1;
        case 2
            xi_x_tArr(2)  = -1;
            eta_x_tArr(2) = -eta_x_tArr(1);
        case 3
            xi_x_tArr(2)  = -eta_x_tArr(1);
            eta_x_tArr(2) =  1;
    end
elseif pointIsInSouth % eta_x_t = -1
    adjacentElements(2) = eNeighbour(e_x,4,1);
    twist = eNeighbour(e_x,4,2);
    switch twist
        case 0
            xi_x_tArr(2)  =  xi_x_tArr(1);
            eta_x_tArr(2) =  1;
        case 1
            xi_x_tArr(2)  =  1;
            eta_x_tArr(2) = -xi_x_tArr(1);
        case 2
            xi_x_tArr(2)  = -xi_x_tArr(1);
            eta_x_tArr(2) = -1;
        case 3
            xi_x_tArr(2)  = -1;
            eta_x_tArr(2) =  xi_x_tArr(1);
    end
end

