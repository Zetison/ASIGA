from rhino3dm import *
from splipy import *
from numpy import *
from splipy.io import g2
import sys

mfile = File3dm.Read(sys.argv[1])

objects = []
for obj in mfile.Objects:
    print(type(obj.Geometry))
    if type(obj.Geometry) is Brep:
        print(len(obj.Geometry.Surfaces))
        print(len(obj.Geometry.Faces))
        for idx in range(0,len(obj.Geometry.Faces)):
            nsrf = obj.Geometry.Faces[idx].UnderlyingSurface().ToNurbsSurface()
            knotsu = [0]
            for i in nsrf.KnotsU:
                knotsu.append(i)
            knotsu.append(knotsu[len(knotsu)-1])
            knotsu[0] = knotsu[1]

            knotsv = [0]
            for i in nsrf.KnotsV:
                knotsv.append(i)
            knotsv.append(knotsv[len(knotsv)-1])
            knotsv[0] = knotsv[1]

            basisu = BSplineBasis(nsrf.OrderU, knotsu, -1)
            basisv = BSplineBasis(nsrf.OrderV, knotsv, -1)
            cpts = []

            cpts = ndarray((nsrf.Points.CountU*nsrf.Points.CountV, 3 + nsrf.IsRational))
            for v in range(0,nsrf.Points.CountV):
                for u in range(0,nsrf.Points.CountU):
                    cpts[u+v*nsrf.Points.CountU,0] = nsrf.Points[u,v].X
                    cpts[u+v*nsrf.Points.CountU,1] = nsrf.Points[u,v].Y
                    cpts[u+v*nsrf.Points.CountU,2] = nsrf.Points[u,v].Z
                    if nsrf.IsRational:
                        cpts[u+v*nsrf.Points.CountU,3] = nsrf.Points[u,v].W

            objects.append(Surface(basisu, basisv, cpts, nsrf.IsRational))

    if type(obj.Geometry) is NurbsCurve:
        knots = [0]
        for i in obj.Geometry.Knots:
            knots.append(i)
        knots[0] = knots[1]
        knots.append(knots[len(knots)-1])
        basis = BSplineBasis(obj.Geometry.Order, knots, -1)
        cpts = []

        cpts = ndarray((len(obj.Geometry.Points), 3 + obj.Geometry.IsRational))
        for u in range(0,len(obj.Geometry.Points)):
            cpts[u,0] = obj.Geometry.Points[u].X
            cpts[u,1] = obj.Geometry.Points[u].Y
            cpts[u,2] = obj.Geometry.Points[u].Z
            if obj.Geometry.IsRational:
                cpts[u,3] = obj.Geometry.Points[u].W

        objects.append(Curve(basis, cpts, obj.Geometry.IsRational))

    if type(obj.Geometry) is NurbsSurface:
        knotsu = [0]
        for i in obj.Geometry.KnotsU:
            knotsu.append(i)
        knotsu.append(knotsu[len(knotsu)-1])
        knotsu[0] = knotsu[1]

        knotsv = [0]
        for i in obj.Geometry.KnotsV:
            knotsv.append(i)
        knotsv.append(knotsv[len(knotsv)-1])
        knotsv[0] = knotsv[1]

        basisu = BSplineBasis(obj.Geometry.OrderU, knotsu, -1)
        basisv = BSplineBasis(obj.Geometry.OrderV, knotsv, -1)
        cpts = []

        nsrf = obj.Geometry

        cpts = ndarray((nsrf.CountU*nsrf.CountV, 3 + nsrf.IsRational))
        for v in range(0,nsrf.Points.CountV):
            for u in range(0,nsrf.Points.CountU):
                cpts[u,v,0] = nsrf.Points[u,v].X
                cpts[u,v,1] = nsrf.Points[u,v].Y
                cpts[u,v,2] = nsrf.Points[u,v].Z
                if nsrf.IsRational:
                    cpts[u,v,3] = nsrf.Points[u,v].W

        objects.append(Surface(basisu, basisv, cpts, obj.Geometry.IsRational))

if len(objects) > 0:
#     print(objects)
    with g2.G2('test.g2') as ofile:
        ofile.write(objects)
