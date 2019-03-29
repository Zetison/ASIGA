function solid = getScordeliLoRoofData(R_i, R_o, L, phi)

solid = getPartedCylinderData(R_i, R_o, L, pi/2-phi, phi);

