function strain = calculateStrainVector(dUdx, dUdy, dUdz)

strain = [dUdx(1);
          dUdy(2);
          dUdz(3);
          dUdz(2) + dUdy(3);
          dUdx(3) + dUdz(1);
          dUdy(1) + dUdx(2)];
          