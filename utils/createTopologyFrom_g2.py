from splipy.io import G2
from splipy import SplineModel

 

model = SplineModel(dimension=3, pardim=2)
with G2('S1.g2') as f:
    model.add(f.read())

 

model.write_ifem('S1_out')
