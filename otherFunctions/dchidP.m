function [dchidtheta, dchidphi] = dchidP(drdX,J,J3)

dchidPvec = drdX*J/J3;

dchidtheta = dchidPvec(1);
dchidphi = dchidPvec(2);      