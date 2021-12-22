function solid = getScordeliLoRoofData(R, R_o, L, phi)
error('Depricated')

solid = getPartedCylinderData(R, R_o, L, pi/2-phi, phi);

