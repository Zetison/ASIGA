function nurbs = getHorseShoeData(R, L, H, D)
error('Depricated')

% 
% Xi = [0 0 0 1 1 2 2 2]/2;
% Eta = [0 0 1 1];
% Zeta = [0 0 0 1 1 2 2 3 3 4 4 4]/4;
% 
% controlPts = zeros(4,5,2,9);
% 
% controlPts(:,1,1,1) = [0                R               0   1                   ];
% controlPts(:,2,1,1) = [-R*(sqrt(2)-1)   R               0   1/2*sqrt(2+sqrt(2))	];
% controlPts(:,3,1,1) = [-R/sqrt(2)       R/sqrt(2)       0   1                   ];
% controlPts(:,4,1,1) = [-R               R*(sqrt(2)-1)   0   1/2*sqrt(2+sqrt(2))	];
% controlPts(:,5,1,1) = [-R               0               0   1                   ];
% controlPts(:,1,2,1) = [0                L               0   1                   ];
% controlPts(:,2,2,1) = [-L/2             L               0   1                   ];
% controlPts(:,3,2,1) = [-L               L               0   1                   ];
% controlPts(:,4,2,1) = [-L               L/2             0   1                   ];
% controlPts(:,5,2,1) = [-L               0               0   1                   ];
% 
% controlPts(:,1,1,2) = [0                R               H/2   1                   ];
% controlPts(:,2,1,2) = [-R*(sqrt(2)-1)   R               H/2   1/2*sqrt(2+sqrt(2)) ];
% controlPts(:,3,1,2) = [-R/sqrt(2)       R/sqrt(2)       H/2   1                   ];
% controlPts(:,4,1,2) = [-R               R*(sqrt(2)-1)   H/2   1/2*sqrt(2+sqrt(2)) ];
% controlPts(:,5,1,2) = [-R               0               H/2   1                   ];
% controlPts(:,1,2,2) = [0                L               H/2   1                   ];
% controlPts(:,2,2,2) = [-L/2             L               H/2   1                   ];
% controlPts(:,3,2,2) = [-L               L               H/2   1                   ];
% controlPts(:,4,2,2) = [-L               L/2             H/2   1                   ];
% controlPts(:,5,2,2) = [-L               0               H/2   1                   ];
% 
% controlPts(:,1,1,3) = [0                R               H   1                   ];
% controlPts(:,2,1,3) = [-R*(sqrt(2)-1)   R               H   1/2*sqrt(2+sqrt(2))	];
% controlPts(:,3,1,3) = [-R/sqrt(2)       R/sqrt(2)       H   1                   ];
% controlPts(:,4,1,3) = [-R               R*(sqrt(2)-1)   H   1/2*sqrt(2+sqrt(2)) ];
% controlPts(:,5,1,3) = [-R               0               H   1                   ];
% controlPts(:,1,2,3) = [0                L               H   1                   ];
% controlPts(:,2,2,3) = [-L/2             L               H   1                   ];
% controlPts(:,3,2,3) = [-L               L               H   1                   ];
% controlPts(:,4,2,3) = [-L               L/2             H   1                   ];
% controlPts(:,5,2,3) = [-L               0               H   1                   ];
% 
% controlPts(:,1,1,4) = [0                R               H+1/sqrt(2)*(R+D)       1/sqrt(2)                   ];
% controlPts(:,2,1,4) = [-R*(sqrt(2)-1)   R               H+1/2*(sqrt(2)*D+2*R)   1/2*sqrt(2+sqrt(2))/sqrt(2)	];
% controlPts(:,3,1,4) = [-R/sqrt(2)       R/sqrt(2)       H+1/2*(sqrt(2)*D+2*R)   1/sqrt(2)                   ];
% controlPts(:,4,1,4) = [-R               R*(sqrt(2)-1)   H+1/2*(sqrt(2)*D+2*R)   1/2*sqrt(2+sqrt(2))/sqrt(2) ];
% controlPts(:,5,1,4) = [-R               0               H+1/sqrt(2)*(R+D)       1/sqrt(2)                   ];
% controlPts(:,1,2,4) = [0                L               H+1/sqrt(2)*(L+D)       1/sqrt(2)                   ];
% controlPts(:,2,2,4) = [-L/2             L               H+1/sqrt(2)*(3*L/2+D)   1/sqrt(2)                   ];
% controlPts(:,3,2,4) = [-L               L               H+1/sqrt(2)*(2*L+D)     1/sqrt(2)                   ];
% controlPts(:,4,2,4) = [-L               L/2             H+1/sqrt(2)*(3*L/2+D)   1/sqrt(2)                   ];
% controlPts(:,5,2,4) = [-L               0               H+1/sqrt(2)*(L+D)       1/sqrt(2)                   ];
% 
% controlPts(:,1,1,5) = [1/2*(R+D)            1/2*(R-D)            	H+1/sqrt(2)*(R+D)     	1                   ];
% controlPts(:,2,1,5) = [D/2-(1/sqrt(2)-1)*R -D/2-(1/sqrt(2)-1)*R    	H+1/2*(sqrt(2)*D+2*R) 	1/2*sqrt(2+sqrt(2))	];
% controlPts(:,3,1,5) = [D/2                  -D/2                    H+1/2*(sqrt(2)*D+2*R)  	1                   ];
% controlPts(:,4,1,5) = [D/2+(1/sqrt(2)-1)*R -D/2+(1/sqrt(2)-1)*R    	H+1/2*(sqrt(2)*D+2*R) 	1/2*sqrt(2+sqrt(2)) ];
% controlPts(:,5,1,5) = [1/2*(D-R)            -1/2*(R+D)              H+1/sqrt(2)*(R+D)    	1                   ];
% controlPts(:,1,2,5) = [(L+D)/2              (L-D)/2                 H+1/sqrt(2)*(L+D)     	1                   ];
% controlPts(:,2,2,5) = [1/4*L+1/2*D          1/4*L-1/2*D             H+1/sqrt(2)*(3*L/2+D)   1                   ];
% controlPts(:,3,2,5) = [D/2                  -D/2                    H+1/sqrt(2)*(2*L+D)     1                   ];
% controlPts(:,4,2,5) = [-L/4+D/2             -L/4-D/2                H+1/sqrt(2)*(3*L/2+D)   1                   ];
% controlPts(:,5,2,5) = [-(L-D)/2             -(L+D)/2                H+1/sqrt(2)*(L+D)       1                   ];
% 
% controlPts(:,1,1,6) = [R+D             	-D                  H+1/sqrt(2)*(R+D)       1/sqrt(2)                   ];
% controlPts(:,2,1,6) = [R+D              -R*(sqrt(2)-1)-D   	H+1/2*(sqrt(2)*D+2*R)   1/2*sqrt(2+sqrt(2))/sqrt(2)	];
% controlPts(:,3,1,6) = [R/sqrt(2)+D    	-R/sqrt(2)-D        H+1/2*(sqrt(2)*D+2*R)   1/sqrt(2)                   ];
% controlPts(:,4,1,6) = [R*(sqrt(2)-1)+D	-R-D                H+1/2*(sqrt(2)*D+2*R)   1/2*sqrt(2+sqrt(2))/sqrt(2) ];
% controlPts(:,5,1,6) = [D                -R-D                H+1/sqrt(2)*(R+D)       1/sqrt(2)                   ];
% controlPts(:,1,2,6) = [L+D             	-D                  H+1/sqrt(2)*(L+D)       1/sqrt(2)                   ];
% controlPts(:,2,2,6) = [L+D              -L/2-D            	H+1/sqrt(2)*(3*L/2+D)   1/sqrt(2)                   ];
% controlPts(:,3,2,6) = [L+D            	-L-D                H+1/sqrt(2)*(2*L+D)     1/sqrt(2)                   ];
% controlPts(:,4,2,6) = [L/2+D           	-L-D                H+1/sqrt(2)*(3*L/2+D)   1/sqrt(2)                   ];
% controlPts(:,5,2,6) = [D                -L-D                H+1/sqrt(2)*(L+D)       1/sqrt(2)                   ];
% 
% controlPts(:,1,1,7) = [R+D              -D                  H   1                   ];
% controlPts(:,2,1,7) = [R+D           	-R*(sqrt(2)-1)-D    H   1/2*sqrt(2+sqrt(2))	];
% controlPts(:,3,1,7) = [R/sqrt(2)+D   	-R/sqrt(2)-D        H   1                   ];
% controlPts(:,4,1,7) = [R*(sqrt(2)-1)+D	-R-D                H   1/2*sqrt(2+sqrt(2))	];
% controlPts(:,5,1,7) = [D                -R-D                H   1                   ];
% controlPts(:,1,2,7) = [L+D              -D                  H   1                   ];
% controlPts(:,2,2,7) = [L+D              -L/2-D          	H   1                   ];
% controlPts(:,3,2,7) = [L+D           	-L-D                H   1                   ];
% controlPts(:,4,2,7) = [L/2+D          	-L-D                H   1                   ];
% controlPts(:,5,2,7) = [D                -L-D                H   1                   ];
% 
% controlPts(:,1,1,8) = [R+D              -D                  H/2   1                   ];
% controlPts(:,2,1,8) = [R+D           	-R*(sqrt(2)-1)-D    H/2   1/2*sqrt(2+sqrt(2)) ];
% controlPts(:,3,1,8) = [R/sqrt(2)+D   	-R/sqrt(2)-D        H/2   1                   ];
% controlPts(:,4,1,8) = [R*(sqrt(2)-1)+D	-R-D                H/2   1/2*sqrt(2+sqrt(2)) ];
% controlPts(:,5,1,8) = [D                -R-D                H/2   1                   ];
% controlPts(:,1,2,8) = [L+D              -D                  H/2   1                   ];
% controlPts(:,2,2,8) = [L+D              -L/2-D          	H/2   1                   ];
% controlPts(:,3,2,8) = [L+D           	-L-D                H/2   1                   ];
% controlPts(:,4,2,8) = [L/2+D          	-L-D                H/2   1                   ];
% controlPts(:,5,2,8) = [D                -L-D                H/2   1                   ];
% 
% controlPts(:,1,1,9) = [R+D              -D                  0   1                   ];
% controlPts(:,2,1,9) = [R+D           	-R*(sqrt(2)-1)-D    0   1/2*sqrt(2+sqrt(2))	];
% controlPts(:,3,1,9) = [R/sqrt(2)+D   	-R/sqrt(2)-D        0   1                   ];
% controlPts(:,4,1,9) = [R*(sqrt(2)-1)+D	-R-D                0   1/2*sqrt(2+sqrt(2))	];
% controlPts(:,5,1,9) = [D                -R-D                0   1                   ];
% controlPts(:,1,2,9) = [L+D              -D                  0   1                   ];
% controlPts(:,2,2,9) = [L+D              -L/2-D          	0   1                   ];
% controlPts(:,3,2,9) = [L+D           	-L-D                0   1                   ];
% controlPts(:,4,2,9) = [L/2+D          	-L-D                0   1                   ];
% controlPts(:,5,2,9) = [D                -L-D                0   1                   ];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% points from Kjetil Johannsen

Xi = [0 0 0 0.5 1 1 1];
Eta = [0 0 0 1 1 1];
Zeta = [0 0 0 1 2 3 3 4 5 6 6 6];

controlPts = zeros(4,4,3,9);

controlPts(:,1,1,1) = [  	-1              0       	0    	1           	];
controlPts(:,2,1,1) = [  	-1              sqrt(2)-1   0    	0.853553391           	];
controlPts(:,3,1,1) = [  	-(sqrt(2)-1)  	1       	0    	0.853553391           	];
controlPts(:,4,1,1) = [  	0               1       	0    	1           	];
controlPts(:,1,2,1) = [  	-2.500      	0       	0    	1           	];
controlPts(:,2,2,1) = [  	-2.437867966   	0.952691193       	0    	1           	];
controlPts(:,3,2,1) = [  	-0.973401872   	2.417157288       	0    	1           	];
controlPts(:,4,2,1) = [  	0               2.500       	0    	1           	];
controlPts(:,1,3,1) = [  	-4              0       	0    	1           	];
controlPts(:,2,3,1) = [  	-4              4       	0    	1           	];
controlPts(:,3,3,1) = [  	-4              4       	0    	1           	];
controlPts(:,4,3,1) = [  	0               4       	0    	1           	];

controlPts(:,1,1,2) = [  	-1              0       	2    	1           	];
controlPts(:,2,1,2) = [  	-1              sqrt(2)-1       	2    	0.853553391           	];
controlPts(:,3,1,2) = [  	-(sqrt(2)-1) 	1       	2    	0.853553391           	];
controlPts(:,4,1,2) = [  	0               1       	2    	1           	];
controlPts(:,1,2,2) = [  	-2.500      	0       	2    	1           	];
controlPts(:,2,2,2) = [  	-2.437867966   	0.952691193       	2    	1           	];
controlPts(:,3,2,2) = [  	-0.973401872  	2.417157288       	2    	1           	];
controlPts(:,4,2,2) = [  	0               2.500       	2    	1           	];
controlPts(:,1,3,2) = [  	-4              0       	2    	1           	];
controlPts(:,2,3,2) = [  	-4              4       	2    	1           	];
controlPts(:,3,3,2) = [  	-4              4       	2    	1           	];
controlPts(:,4,3,2) = [  	0               4       	2    	1           	];

controlPts(:,1,1,3) = [  	-1              0       	4    	1           	];
controlPts(:,2,1,3) = [  	-1              sqrt(2)-1       	4    	0.853553391           	];
controlPts(:,3,1,3) = [  	-(sqrt(2)-1)   	1       	4    	0.853553391           	];
controlPts(:,4,1,3) = [  	0               1     	4    	1           	];
controlPts(:,1,2,3) = [  	-2.500      	0       	4    	1           	];
controlPts(:,2,2,3) = [  	-2.437867966  	0.952691193       	4    	1           	];
controlPts(:,3,2,3) = [  	-0.973401872   	2.417157288       	4    	1           	];
controlPts(:,4,2,3) = [  	0               2.500       	4    	1           	];
controlPts(:,1,3,3) = [  	-4              0       	4    	1           	];
controlPts(:,2,3,3) = [  	-4              4       	4    	1           	];
controlPts(:,3,3,3) = [  	-4              4       	4    	1           	];
controlPts(:,4,3,3) = [  	0               4       	4    	1           	];

controlPts(:,1,1,4) = [  	-1              0       	7.0178    	1/sqrt(2)           	];
controlPts(:,2,1,4) = [  	-1              sqrt(2)-1       	7.4393    	0.603553391           	];
controlPts(:,3,1,4) = [  	-(sqrt(2)-1)   	1       	7.4393    	0.603553391           	];
controlPts(:,4,1,4) = [  	0               1       	7.0178    	1/sqrt(2)           	];
controlPts(:,1,2,4) = [  	-2.500      	0       	8.5444    	1/sqrt(2)           	];
controlPts(:,2,2,4) = [  	-2.4378679663  	0.952691193       	9.450800001    	1/sqrt(2)           	];
controlPts(:,3,2,4) = [  	-0.973401872  	2.417157288       	9.450800001    	1/sqrt(2)           	];
controlPts(:,4,2,4) = [  	0               2.500       	8.5444    	1/sqrt(2)           	];
controlPts(:,1,3,4) = [  	-4              0       	10.0711    	1/sqrt(2)           	];
controlPts(:,2,3,4) = [  	-4              4       	14.4142    	1/sqrt(2)           	];
controlPts(:,3,3,4) = [  	-4              4       	14.4142    	1/sqrt(2)           	];
controlPts(:,4,3,4) = [  	0               4       	10.0711    	1/sqrt(2)           	];

controlPts(:,1,1,5) = [  	0.2071      	-1.2071       	7.0178    	1           	];
controlPts(:,2,1,5) = [  	sqrt(2)-1      	-1       	7.4393    	0.853553391           	];
controlPts(:,3,1,5) = [  	1      	-(sqrt(2)-1)       	7.4393    	0.853553391           	];
controlPts(:,4,1,5) = [  	1.2071      	-0.2071       	7.0178    	1           	];
controlPts(:,1,2,5) = [  	-0.542893219      	-1.957106781       	8.5444    	1           	];
controlPts(:,2,2,5) = [  	-0.035481605      	-1.449695168       	9.4508    	1           	];
controlPts(:,3,2,5) = [  	1.428984489      	0.014770927       	9.4508    	1           	];
controlPts(:,4,2,5) = [  	1.957106781      	0.542893219       	8.5444    	1           	];
controlPts(:,1,3,5) = [  	-1.292893219      	-2.707106781       	10.0711    	1           	];
controlPts(:,2,3,5) = [  	1/sqrt(2)      	-1/sqrt(2)       	14.4142    	1           	];
controlPts(:,3,3,5) = [  	1/sqrt(2)      	-1/sqrt(2)       	14.4142    	1           	];
controlPts(:,4,3,5) = [  	2.707106781      	1.292893219       	10.0711    	1           	];

controlPts(:,1,1,6) = [  	1      	-2       	7.0178    	1/sqrt(2)           	];
controlPts(:,2,1,6) = [  	sqrt(2)      	-2       	7.4393    	0.603553391           	];
controlPts(:,3,1,6) = [  	2      	-sqrt(2)       	7.4393    	0.603553391           	];
controlPts(:,4,1,6) = [  	2      	-1       	7.0178    	1/sqrt(2)           	];
controlPts(:,1,2,6) = [  	1      	-3.5       	8.5444    	1/sqrt(2)           	];
controlPts(:,2,2,6) = [  	1.952691193      	-3.4378679663       	9.450800001    	1/sqrt(2)           	];
controlPts(:,3,2,6) = [  	3.417157288      	-1.973401871999994       	9.450800001    	1/sqrt(2)           	];
controlPts(:,4,2,6) = [  	3.5      	-1       	8.5444    	1/sqrt(2)           	];
controlPts(:,1,3,6) = [  	1      	-5       	10.0711    	1/sqrt(2)           	];
controlPts(:,2,3,6) = [  	5      	-5       	14.4142    	1/sqrt(2)           	];
controlPts(:,3,3,6) = [  	5      	-5       	14.4142    	1/sqrt(2)           	];
controlPts(:,4,3,6) = [  	5      	-1       	10.0711    	1/sqrt(2)           	];

controlPts(:,1,1,7) = [  	1      	-2       	4    	1           	];
controlPts(:,2,1,7) = [  	sqrt(2)      	-2       	4    	0.853553391           	];
controlPts(:,3,1,7) = [  	2      	-sqrt(2)       	4    	0.853553391           	];
controlPts(:,4,1,7) = [  	2      	-1       	4    	1           	];
controlPts(:,1,2,7) = [  	1      	-3.5       	4    	1           	];
controlPts(:,2,2,7) = [  	1.952691193      	-3.437867966       	4    	1           	];
controlPts(:,3,2,7) = [  	3.417157288      	-1.973401872       	4    	1           	];
controlPts(:,4,2,7) = [  	3.5      	-1       	4    	1           	];
controlPts(:,1,3,7) = [  	1      	-5       	4    	1           	];
controlPts(:,2,3,7) = [  	5      	-5       	4    	1           	];
controlPts(:,3,3,7) = [  	5      	-5       	4    	1           	];
controlPts(:,4,3,7) = [  	5      	-1       	4    	1           	];

controlPts(:,1,1,8) = [  	1      	-2       	2    	1           	];
controlPts(:,2,1,8) = [  	sqrt(2)      	-2       	2    	0.853553391           	];
controlPts(:,3,1,8) = [  	2      	-sqrt(2)       	2    	0.853553391           	];
controlPts(:,4,1,8) = [  	2      	-1       	2    	1           	];
controlPts(:,1,2,8) = [  	1      	-3.5       	2    	1           	];
controlPts(:,2,2,8) = [  	1.952691193      	-3.437867966       	2    	1           	];
controlPts(:,3,2,8) = [  	3.417157288      	-1.973401872       	2    	1           	];
controlPts(:,4,2,8) = [  	3.5      	-1       	2    	1           	];
controlPts(:,1,3,8) = [  	1      	-5       	2    	1           	];
controlPts(:,2,3,8) = [  	5      	-5       	2    	1           	];
controlPts(:,3,3,8) = [  	5      	-5       	2    	1           	];
controlPts(:,4,3,8) = [  	5      	-1       	2    	1           	];

controlPts(:,1,1,9) = [  	1      	-2       	0    	1           	];
controlPts(:,2,1,9) = [  	sqrt(2)      	-2       	0    	0.853553391           	];
controlPts(:,3,1,9) = [  	2      	-sqrt(2)       	0    	0.853553391           	];
controlPts(:,4,1,9) = [  	2      	-1       	0    	1           	];
controlPts(:,1,2,9) = [  	1      	-3.5       	0    	1           	];
controlPts(:,2,2,9) = [  	1.952691193      	-3.437867966       	0    	1           	];
controlPts(:,3,2,9) = [  	3.417157288      	-1.973401872       	0    	1           	];
controlPts(:,4,2,9) = [  	3.5      	-1       	0    	1           	];
controlPts(:,1,3,9) = [  	1      	-5       	0    	1           	];
controlPts(:,2,3,9) = [  	5      	-5       	0    	1           	];
controlPts(:,3,3,9) = [  	5      	-5       	0    	1           	];
controlPts(:,4,3,9) = [  	5      	-1       	0    	1           	];

controlPts(3,:,:,:) = controlPts(3,:,:,:)*-1;

nurbs = createNURBSobject(controlPts,{Xi, Eta, Zeta});