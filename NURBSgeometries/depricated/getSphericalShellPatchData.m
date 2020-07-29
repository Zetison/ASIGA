function nurbs = getSphericalShellPatchData(R_o,t)
error('Use getEllipsoidData()')
% data from Cobb1988tts
Xi = [0 0 0 0 0 1 1 1 1 1];
Eta = [0 0 0 0 0 1 1 1 1 1];
Zeta = [0 0 1 1];
R_i = R_o-t;

nurbs = cell(1,6);
if isa(R_o,'sym')
    sr2 = sqrt(sym('2'));                                                        
    sr3 = sqrt(sym('3'));                                                        
    sr6 = sqrt(sym('6'));   
else
    sr2 = sqrt(2);                                                        
    sr3 = sqrt(3);                                                        
    sr6 = sqrt(6);   
end
controlPts = zeros(4,5,5,2,class(R_o)); 

if true
    controlPts(:,1,1,1) = [4*(1-sr3)    4*(1-sr3)       4*(sr3-1)       4*(3-sr3)];
    controlPts(:,2,1,1) = [-sr2         sr2*(sr3-4)     sr2*(4-sr3)     sr2*(3*sr3-2)];
    controlPts(:,3,1,1) = [0            4*(1-2*sr3)/3   4*(2*sr3-1)/3	4*(5-sr3)/3]; 
    controlPts(:,2,2,1) = [-(3*sr3-2)/2	(2-3*sr3)/2     (sr3+6)/2       (sr3+6)/2];
    controlPts(:,3,2,1) = [0            sr2*(2*sr3-7)/3	5*sr6/3         sr2*(sr3+6)/3];
    controlPts(:,3,3,1) = [0            0               4*(5-sr3)/3     4*(5*sr3-1)/9];
    pairs = [1,2; 1,3; 2,3];
    for l = 1:3
        i = pairs(l,1);
        j = pairs(l,2);
        controlPts(1,i,j,1) = controlPts(2,j,i,1);
        controlPts(2,i,j,1) = controlPts(1,j,i,1);
        controlPts(3:4,i,j,1) = controlPts(3:4,j,i,1);
    end
    for i = 4:5
        for j = 1:3
            controlPts(1,i,j,1) = -controlPts(1,6-i,j,1);
            controlPts(2:4,i,j,1) = controlPts(2:4,6-i,j,1);
        end
    end
    for i = 1:5
        for j = 4:5
            controlPts(2,i,j,1) = -controlPts(2,i,6-j,1);
            controlPts([1,3,4],i,j,1) = controlPts([1,3,4],i,6-j,1);
        end
    end
else              
    controlPts(:,1,1,1) = [4*(1-sr3)    4*(1-sr3)       4*(1-sr3)       4*(3-sr3)];
    controlPts(:,2,1,1) = [-sr2         sr2*(sr3-4)     sr2*(sr3-4)     sr2*(3*sr3-2)];
    controlPts(:,3,1,1) = [0            4*(1-2*sr3)/3   4*(1-2*sr3)/3	4*(5-sr3)/3];
    controlPts(:,4,1,1) = [sr2          sr2*(sr3-4)     sr2*(sr3-4)     sr2*(3*sr3-2)];
    controlPts(:,5,1,1) = [4*(sr3-1)	4*(1-sr3)       4*(1-sr3)       4*(3-sr3)];

    controlPts(:,1,2,1) = [-sr2*(4-sr3)	-sr2            sr2*(sr3-4) sr2*(3*sr3-2)];
    controlPts(:,2,2,1) = [-(3*sr3-2)/2	(2-3*sr3)/2     -(sr3+6)/2  (sr3+6)/2];
    controlPts(:,3,2,1) = [0            sr2*(2*sr3-7)/3	-5*sr6/3    sr2*(sr3+6)/3];
    controlPts(:,4,2,1) = [(3*sr3-2)/2	(2-3*sr3)/2     -(sr3+6)/2  (sr3+6)/2];
    controlPts(:,5,2,1) = [sr2*(4-sr3)	-sr2            sr2*(sr3-4) sr2*(3*sr3-2)];

    controlPts(:,1,3,1) = [-4*(2*sr3-1)/3	0           4*(1-2*sr3)/3   4*(5-sr3)/3];
    controlPts(:,2,3,1) = [-sr2*(7-2*sr3)/3	0           -5*sr6/3        sr2*(sr3+6)/3];
    controlPts(:,3,3,1) = [0                0           4*(sr3-5)/3     4*(5*sr3-1)/9];
    controlPts(:,4,3,1) = [sr2*(7-2*sr3)/3	0           -5*sr6/3        sr2*(sr3+6)/3];
    controlPts(:,5,3,1) = [4*(2*sr3-1)/3  	0           4*(1-2*sr3)/3   4*(5-sr3)/3];

    controlPts(:,1,4,1) = [-sr2*(4-sr3) 	sr2                 sr2*(sr3-4)	sr2*(3*sr3-2)];
    controlPts(:,2,4,1) = [-(3*sr3-2)/2     -(2-3*sr3)/2        -(sr3+6)/2	(sr3+6)/2];
    controlPts(:,3,4,1) = [0                -sr2*(2*sr3-7)/3 	-5*sr6/3    sr2*(sr3+6)/3];
    controlPts(:,4,4,1) = [(3*sr3-2)/2      -(2-3*sr3)/2        -(sr3+6)/2  (sr3+6)/2];
    controlPts(:,5,4,1) = [sr2*(4-sr3)  	sr2                 sr2*(sr3-4) sr2*(3*sr3-2)];

    controlPts(:,1,5,1) = [4*(1-sr3) 	-4*(1-sr3)      4*(1-sr3)       4*(3-sr3)];
    controlPts(:,2,5,1) = [-sr2         -sr2*(sr3-4)    sr2*(sr3-4)     sr2*(3*sr3-2)];
    controlPts(:,3,5,1) = [0            -4*(1-2*sr3)/3  4*(1-2*sr3)/3	4*(5-sr3)/3];
    controlPts(:,4,5,1) = [sr2          -sr2*(sr3-4)    sr2*(sr3-4)     sr2*(3*sr3-2)];
    controlPts(:,5,5,1) = [4*(sr3-1)  	-4*(1-sr3)      4*(1-sr3)       4*(3-sr3)];
    controlPts(3,:,:,:) = -controlPts(3,:,:,:); % fix the orientation such that the normal vector points outwards
end


for i = 1:5
    for j = 1:5
        controlPts(1:3,i,j,1) = controlPts(1:3,i,j,1)/controlPts(4,i,j,1);
    end
end  
controlPts(1:3,:,:,2) = controlPts(1:3,:,:,1)*R_o;  
controlPts(4,:,:,2) = controlPts(4,:,:,1);  
controlPts(1:3,:,:,1) = controlPts(1:3,:,:,1)*R_i; 
nurbs{1} = createNURBSobject(controlPts,{Xi, Eta, Zeta}); 
