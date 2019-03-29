function nurbs = getInfinitePlateWithCircularHoleData(R, L, meshtype)

switch meshtype
    case 1
        Xi = [0 0 0 0.5 1 1 1];
        Eta = [0 0 0 1 1 1];
        H = (R+L)/2;
        controlPts = zeros(3,4,3);

        controlPts(:,1,1) = [-R             0               1           	];
        controlPts(:,2,1) = [-R             R*(sqrt(2)-1)	(1+1/sqrt(2))/2 ];
        controlPts(:,3,1) = [ R*(1-sqrt(2)) R               (1+1/sqrt(2))/2 ];
        controlPts(:,4,1) = [ 0            	R               1           	];

        controlPts(:,1,2) = [-H             0               1           	];
        controlPts(:,2,2) = [-H             0.75*R         	1           	];
        controlPts(:,3,2) = [-0.75*R        H           	1           	];
        controlPts(:,4,2) = [ 0             H       		1           	];

        controlPts(:,1,3) = [-L              0       	   	1           	];
        controlPts(:,2,3) = [-L              L           	1           	];
        controlPts(:,3,3) = [-L              L       	   	1           	];
        controlPts(:,4,3) = [ 0              L           	1           	];


        nurbs = createNURBSobject(controlPts,{Xi, Eta});
    case 2
        Xi = [0 0 0 0.5 0.5 1 1 1];
        Eta = [0 0 1 1];
        
        controlPts = zeros(3,5,2);

        controlPts(:,1,1) = [-R               0           	1                   ];
        controlPts(:,2,1) = [-R               R*(sqrt(2)-1)	sqrt(2+sqrt(2))/2	];
        controlPts(:,3,1) = [-R/sqrt(2)       R/sqrt(2)   	1                   ];
        controlPts(:,4,1) = [-R*(sqrt(2)-1)   R          	sqrt(2+sqrt(2))/2	];
        controlPts(:,5,1) = [0                R             1                   ];
        
        controlPts(:,1,2) = [-L               0           	1                   ];
        controlPts(:,2,2) = [-L               L/2         	1                   ];
        controlPts(:,3,2) = [-L               L          	1                   ];
        controlPts(:,4,2) = [-L/2             L           	1                   ];
        controlPts(:,5,2) = [0                L             1                   ];
        
        
        nurbs = createNURBSobject(controlPts,{Xi, Eta});
        nurbs = elevateNURBSdegree(nurbs,[0 1]);
end
