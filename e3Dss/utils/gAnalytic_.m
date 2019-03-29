function dp = gAnalytic_(v,options)

data = e3Dss(v,options);
dp = [data(1).dpdx, data(1).dpdy, data(1).dpdz];

