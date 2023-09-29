# generated using paraview version 5.8.1-RC2-1072-g4e1578c1bd
#
# To ensure correct image size when batch processing, please search
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# get animation scene
animationScene1 = GetAnimationScene()

# get the time-keeper
timeKeeper1 = GetTimeKeeper()

# create a new 'XML Unstructured Grid Reader'
s1_M5_IE_SHBC_1vtu = XMLUnstructuredGridReader(registrationName='S1_M5_IE_SHBC_1.vtu', FileName=['/home/zetison/results/ASIGA/S1/paraviewResults/S1_M5_IE_SHBC_1.vtu'])
s1_M5_IE_SHBC_1vtu.PointArrayStatus = ['Displacement', 'Scalar field', 'Total scalar field (real)', 'Total scalar field (abs)', 'P_inc', 'Analytic', 'Error', 'Error Gradient', 'Error in Energy']

# create a new 'XML Unstructured Grid Reader'
s1_M5_IE_SHBC_1_meshvtu = XMLUnstructuredGridReader(registrationName='S1_M5_IE_SHBC_1_mesh.vtu', FileName=['/home/zetison/results/ASIGA/S1/paraviewResults/S1_M5_IE_SHBC_1_mesh.vtu'])
s1_M5_IE_SHBC_1_meshvtu.PointArrayStatus = ['Displacement']

# create a new 'XML Unstructured Grid Reader'
s1_M5_PML_SHBC_1vtu = XMLUnstructuredGridReader(registrationName='S1_M5_PML_SHBC_1.vtu', FileName=['/home/zetison/results/ASIGA/S1/paraviewResults/S1_M5_PML_SHBC_1.vtu'])
s1_M5_PML_SHBC_1vtu.PointArrayStatus = ['Displacement', 'Scalar field', 'Total scalar field (real)', 'Total scalar field (abs)', 'P_inc', 'Analytic', 'Error', 'Error Gradient', 'Error in Energy']

# create a new 'XML Unstructured Grid Reader'
s1_M5_PML_SHBC_1_meshvtu = XMLUnstructuredGridReader(registrationName='S1_M5_PML_SHBC_1_mesh.vtu', FileName=['/home/zetison/results/ASIGA/S1/paraviewResults/S1_M5_PML_SHBC_1_mesh.vtu'])
s1_M5_PML_SHBC_1_meshvtu.PointArrayStatus = ['Displacement']

# Properties modified on s1_M5_IE_SHBC_1_meshvtu
s1_M5_IE_SHBC_1_meshvtu.TimeArray = 'None'

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1450, 867]

# get layout
layout1 = GetLayout()

# show data in view
s1_M5_IE_SHBC_1_meshvtuDisplay = Show(s1_M5_IE_SHBC_1_meshvtu, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
s1_M5_IE_SHBC_1_meshvtuDisplay.Representation = 'Surface'
s1_M5_IE_SHBC_1_meshvtuDisplay.ColorArrayName = [None, '']
s1_M5_IE_SHBC_1_meshvtuDisplay.SelectTCoordArray = 'None'
s1_M5_IE_SHBC_1_meshvtuDisplay.SelectNormalArray = 'None'
s1_M5_IE_SHBC_1_meshvtuDisplay.SelectTangentArray = 'None'
s1_M5_IE_SHBC_1_meshvtuDisplay.OSPRayScaleArray = 'Displacement'
s1_M5_IE_SHBC_1_meshvtuDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
s1_M5_IE_SHBC_1_meshvtuDisplay.SelectOrientationVectors = 'Displacement'
s1_M5_IE_SHBC_1_meshvtuDisplay.ScaleFactor = 0.30000000000000004
s1_M5_IE_SHBC_1_meshvtuDisplay.SelectScaleArray = 'Displacement'
s1_M5_IE_SHBC_1_meshvtuDisplay.GlyphType = 'Arrow'
s1_M5_IE_SHBC_1_meshvtuDisplay.GlyphTableIndexArray = 'Displacement'
s1_M5_IE_SHBC_1_meshvtuDisplay.GaussianRadius = 0.015
s1_M5_IE_SHBC_1_meshvtuDisplay.SetScaleArray = ['POINTS', 'Displacement']
s1_M5_IE_SHBC_1_meshvtuDisplay.ScaleTransferFunction = 'PiecewiseFunction'
s1_M5_IE_SHBC_1_meshvtuDisplay.OpacityArray = ['POINTS', 'Displacement']
s1_M5_IE_SHBC_1_meshvtuDisplay.OpacityTransferFunction = 'PiecewiseFunction'
s1_M5_IE_SHBC_1_meshvtuDisplay.DataAxesGrid = 'GridAxesRepresentation'
s1_M5_IE_SHBC_1_meshvtuDisplay.PolarAxes = 'PolarAxesRepresentation'
s1_M5_IE_SHBC_1_meshvtuDisplay.ScalarOpacityUnitDistance = 0.3121966330432184
s1_M5_IE_SHBC_1_meshvtuDisplay.OpacityArrayName = ['POINTS', 'Displacement']
s1_M5_IE_SHBC_1_meshvtuDisplay.SelectInputVectors = ['POINTS', 'Displacement']
s1_M5_IE_SHBC_1_meshvtuDisplay.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
s1_M5_IE_SHBC_1_meshvtuDisplay.OSPRayScaleFunction.Points = [-0.0155398497892645, 0.0, 0.5, 0.0, -0.0009743210460688503, 0.5, 0.5, 0.0, 0.0135912076971268, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
s1_M5_IE_SHBC_1_meshvtuDisplay.ScaleTransferFunction.Points = [-4.13450632764285e-11, 0.0, 0.5, 0.0, 1.43459622751001e-14, 0.5, 0.5, 0.0, 4.13737552009787e-11, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
s1_M5_IE_SHBC_1_meshvtuDisplay.OpacityTransferFunction.Points = [-4.13450632764285e-11, 0.0, 0.5, 0.0, 1.43459622751001e-14, 0.5, 0.5, 0.0, 4.13737552009787e-11, 1.0, 0.5, 0.0]

# reset view to fit data
renderView1.ResetCamera()

# Properties modified on s1_M5_IE_SHBC_1vtu
s1_M5_IE_SHBC_1vtu.TimeArray = 'None'

# show data in view
s1_M5_IE_SHBC_1vtuDisplay = Show(s1_M5_IE_SHBC_1vtu, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
s1_M5_IE_SHBC_1vtuDisplay.Representation = 'Surface'
s1_M5_IE_SHBC_1vtuDisplay.ColorArrayName = [None, '']
s1_M5_IE_SHBC_1vtuDisplay.SelectTCoordArray = 'None'
s1_M5_IE_SHBC_1vtuDisplay.SelectNormalArray = 'None'
s1_M5_IE_SHBC_1vtuDisplay.SelectTangentArray = 'None'
s1_M5_IE_SHBC_1vtuDisplay.OSPRayScaleArray = 'Analytic'
s1_M5_IE_SHBC_1vtuDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
s1_M5_IE_SHBC_1vtuDisplay.SelectOrientationVectors = 'Analytic'
s1_M5_IE_SHBC_1vtuDisplay.ScaleFactor = 0.30000000000000004
s1_M5_IE_SHBC_1vtuDisplay.SelectScaleArray = 'Analytic'
s1_M5_IE_SHBC_1vtuDisplay.GlyphType = 'Arrow'
s1_M5_IE_SHBC_1vtuDisplay.GlyphTableIndexArray = 'Analytic'
s1_M5_IE_SHBC_1vtuDisplay.GaussianRadius = 0.015
s1_M5_IE_SHBC_1vtuDisplay.SetScaleArray = ['POINTS', 'Analytic']
s1_M5_IE_SHBC_1vtuDisplay.ScaleTransferFunction = 'PiecewiseFunction'
s1_M5_IE_SHBC_1vtuDisplay.OpacityArray = ['POINTS', 'Analytic']
s1_M5_IE_SHBC_1vtuDisplay.OpacityTransferFunction = 'PiecewiseFunction'
s1_M5_IE_SHBC_1vtuDisplay.DataAxesGrid = 'GridAxesRepresentation'
s1_M5_IE_SHBC_1vtuDisplay.PolarAxes = 'PolarAxesRepresentation'
s1_M5_IE_SHBC_1vtuDisplay.ScalarOpacityUnitDistance = 0.1805890446558368
s1_M5_IE_SHBC_1vtuDisplay.OpacityArrayName = ['POINTS', 'Analytic']
s1_M5_IE_SHBC_1vtuDisplay.SelectInputVectors = ['POINTS', 'Displacement']
s1_M5_IE_SHBC_1vtuDisplay.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
s1_M5_IE_SHBC_1vtuDisplay.OSPRayScaleFunction.Points = [-0.0155398497892645, 0.0, 0.5, 0.0, -0.0009743210460688503, 0.5, 0.5, 0.0, 0.0135912076971268, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
s1_M5_IE_SHBC_1vtuDisplay.ScaleTransferFunction.Points = [-1.32295140808951, 0.0, 0.5, 0.0, 0.3214089293358051, 0.5, 0.5, 0.0, 1.9657692667611202, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
s1_M5_IE_SHBC_1vtuDisplay.OpacityTransferFunction.Points = [-1.32295140808951, 0.0, 0.5, 0.0, 0.3214089293358051, 0.5, 0.5, 0.0, 1.9657692667611202, 1.0, 0.5, 0.0]

# Properties modified on s1_M5_PML_SHBC_1_meshvtu
s1_M5_PML_SHBC_1_meshvtu.TimeArray = 'None'

# show data in view
s1_M5_PML_SHBC_1_meshvtuDisplay = Show(s1_M5_PML_SHBC_1_meshvtu, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
s1_M5_PML_SHBC_1_meshvtuDisplay.Representation = 'Surface'
s1_M5_PML_SHBC_1_meshvtuDisplay.ColorArrayName = [None, '']
s1_M5_PML_SHBC_1_meshvtuDisplay.SelectTCoordArray = 'None'
s1_M5_PML_SHBC_1_meshvtuDisplay.SelectNormalArray = 'None'
s1_M5_PML_SHBC_1_meshvtuDisplay.SelectTangentArray = 'None'
s1_M5_PML_SHBC_1_meshvtuDisplay.OSPRayScaleArray = 'Displacement'
s1_M5_PML_SHBC_1_meshvtuDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
s1_M5_PML_SHBC_1_meshvtuDisplay.SelectOrientationVectors = 'Displacement'
s1_M5_PML_SHBC_1_meshvtuDisplay.ScaleFactor = 0.30000000000000004
s1_M5_PML_SHBC_1_meshvtuDisplay.SelectScaleArray = 'Displacement'
s1_M5_PML_SHBC_1_meshvtuDisplay.GlyphType = 'Arrow'
s1_M5_PML_SHBC_1_meshvtuDisplay.GlyphTableIndexArray = 'Displacement'
s1_M5_PML_SHBC_1_meshvtuDisplay.GaussianRadius = 0.015
s1_M5_PML_SHBC_1_meshvtuDisplay.SetScaleArray = ['POINTS', 'Displacement']
s1_M5_PML_SHBC_1_meshvtuDisplay.ScaleTransferFunction = 'PiecewiseFunction'
s1_M5_PML_SHBC_1_meshvtuDisplay.OpacityArray = ['POINTS', 'Displacement']
s1_M5_PML_SHBC_1_meshvtuDisplay.OpacityTransferFunction = 'PiecewiseFunction'
s1_M5_PML_SHBC_1_meshvtuDisplay.DataAxesGrid = 'GridAxesRepresentation'
s1_M5_PML_SHBC_1_meshvtuDisplay.PolarAxes = 'PolarAxesRepresentation'
s1_M5_PML_SHBC_1_meshvtuDisplay.ScalarOpacityUnitDistance = 0.2228352968694388
s1_M5_PML_SHBC_1_meshvtuDisplay.OpacityArrayName = ['POINTS', 'Displacement']
s1_M5_PML_SHBC_1_meshvtuDisplay.SelectInputVectors = ['POINTS', 'Displacement']
s1_M5_PML_SHBC_1_meshvtuDisplay.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
s1_M5_PML_SHBC_1_meshvtuDisplay.OSPRayScaleFunction.Points = [-0.0155398497892645, 0.0, 0.5, 0.0, -0.0009743210460688503, 0.5, 0.5, 0.0, 0.0135912076971268, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
s1_M5_PML_SHBC_1_meshvtuDisplay.ScaleTransferFunction.Points = [-4.08875818923397e-11, 0.0, 0.5, 0.0, 1.4187224806502317e-14, 0.5, 0.5, 0.0, 4.09159563419527e-11, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
s1_M5_PML_SHBC_1_meshvtuDisplay.OpacityTransferFunction.Points = [-4.08875818923397e-11, 0.0, 0.5, 0.0, 1.4187224806502317e-14, 0.5, 0.5, 0.0, 4.09159563419527e-11, 1.0, 0.5, 0.0]

# Properties modified on s1_M5_PML_SHBC_1vtu
s1_M5_PML_SHBC_1vtu.TimeArray = 'None'

# show data in view
s1_M5_PML_SHBC_1vtuDisplay = Show(s1_M5_PML_SHBC_1vtu, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
s1_M5_PML_SHBC_1vtuDisplay.Representation = 'Surface'
s1_M5_PML_SHBC_1vtuDisplay.ColorArrayName = [None, '']
s1_M5_PML_SHBC_1vtuDisplay.SelectTCoordArray = 'None'
s1_M5_PML_SHBC_1vtuDisplay.SelectNormalArray = 'None'
s1_M5_PML_SHBC_1vtuDisplay.SelectTangentArray = 'None'
s1_M5_PML_SHBC_1vtuDisplay.OSPRayScaleArray = 'Analytic'
s1_M5_PML_SHBC_1vtuDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
s1_M5_PML_SHBC_1vtuDisplay.SelectOrientationVectors = 'Analytic'
s1_M5_PML_SHBC_1vtuDisplay.ScaleFactor = 0.30000000000000004
s1_M5_PML_SHBC_1vtuDisplay.SelectScaleArray = 'Analytic'
s1_M5_PML_SHBC_1vtuDisplay.GlyphType = 'Arrow'
s1_M5_PML_SHBC_1vtuDisplay.GlyphTableIndexArray = 'Analytic'
s1_M5_PML_SHBC_1vtuDisplay.GaussianRadius = 0.015
s1_M5_PML_SHBC_1vtuDisplay.SetScaleArray = ['POINTS', 'Analytic']
s1_M5_PML_SHBC_1vtuDisplay.ScaleTransferFunction = 'PiecewiseFunction'
s1_M5_PML_SHBC_1vtuDisplay.OpacityArray = ['POINTS', 'Analytic']
s1_M5_PML_SHBC_1vtuDisplay.OpacityTransferFunction = 'PiecewiseFunction'
s1_M5_PML_SHBC_1vtuDisplay.DataAxesGrid = 'GridAxesRepresentation'
s1_M5_PML_SHBC_1vtuDisplay.PolarAxes = 'PolarAxesRepresentation'
s1_M5_PML_SHBC_1vtuDisplay.ScalarOpacityUnitDistance = 0.12889829395335262
s1_M5_PML_SHBC_1vtuDisplay.OpacityArrayName = ['POINTS', 'Analytic']
s1_M5_PML_SHBC_1vtuDisplay.SelectInputVectors = ['POINTS', 'Displacement']
s1_M5_PML_SHBC_1vtuDisplay.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
s1_M5_PML_SHBC_1vtuDisplay.OSPRayScaleFunction.Points = [-0.0155398497892645, 0.0, 0.5, 0.0, -0.0009743210460688503, 0.5, 0.5, 0.0, 0.0135912076971268, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
s1_M5_PML_SHBC_1vtuDisplay.ScaleTransferFunction.Points = [-1.40464062906829, 0.0, 0.5, 0.0, 0.2805643188464151, 0.5, 0.5, 0.0, 1.9657692667611202, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
s1_M5_PML_SHBC_1vtuDisplay.OpacityTransferFunction.Points = [-1.40464062906829, 0.0, 0.5, 0.0, 0.2805643188464151, 0.5, 0.5, 0.0, 1.9657692667611202, 1.0, 0.5, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# split cell
layout1.SplitHorizontal(0, 0.5)

# set active view
SetActiveView(None)

# Create a new 'Render View'
renderView2 = CreateView('RenderView')
renderView2.AxesGrid = 'GridAxes3DActor'
renderView2.StereoType = 'Crystal Eyes'
renderView2.CameraFocalDisk = 1.0
# uncomment following to set a specific view size
# renderView2.ViewSize = [400, 400]

# assign view to a particular cell in the layout
AssignViewToLayout(view=renderView2, layout=layout1, hint=2)

# set active source
SetActiveSource(s1_M5_PML_SHBC_1vtu)

# show data in view
s1_M5_PML_SHBC_1vtuDisplay_1 = Show(s1_M5_PML_SHBC_1vtu, renderView2, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
s1_M5_PML_SHBC_1vtuDisplay_1.Representation = 'Surface'
s1_M5_PML_SHBC_1vtuDisplay_1.ColorArrayName = [None, '']
s1_M5_PML_SHBC_1vtuDisplay_1.SelectTCoordArray = 'None'
s1_M5_PML_SHBC_1vtuDisplay_1.SelectNormalArray = 'None'
s1_M5_PML_SHBC_1vtuDisplay_1.SelectTangentArray = 'None'
s1_M5_PML_SHBC_1vtuDisplay_1.OSPRayScaleArray = 'Analytic'
s1_M5_PML_SHBC_1vtuDisplay_1.OSPRayScaleFunction = 'PiecewiseFunction'
s1_M5_PML_SHBC_1vtuDisplay_1.SelectOrientationVectors = 'Analytic'
s1_M5_PML_SHBC_1vtuDisplay_1.ScaleFactor = 0.30000000000000004
s1_M5_PML_SHBC_1vtuDisplay_1.SelectScaleArray = 'Analytic'
s1_M5_PML_SHBC_1vtuDisplay_1.GlyphType = 'Arrow'
s1_M5_PML_SHBC_1vtuDisplay_1.GlyphTableIndexArray = 'Analytic'
s1_M5_PML_SHBC_1vtuDisplay_1.GaussianRadius = 0.015
s1_M5_PML_SHBC_1vtuDisplay_1.SetScaleArray = ['POINTS', 'Analytic']
s1_M5_PML_SHBC_1vtuDisplay_1.ScaleTransferFunction = 'PiecewiseFunction'
s1_M5_PML_SHBC_1vtuDisplay_1.OpacityArray = ['POINTS', 'Analytic']
s1_M5_PML_SHBC_1vtuDisplay_1.OpacityTransferFunction = 'PiecewiseFunction'
s1_M5_PML_SHBC_1vtuDisplay_1.DataAxesGrid = 'GridAxesRepresentation'
s1_M5_PML_SHBC_1vtuDisplay_1.PolarAxes = 'PolarAxesRepresentation'
s1_M5_PML_SHBC_1vtuDisplay_1.ScalarOpacityUnitDistance = 0.12889829395335262
s1_M5_PML_SHBC_1vtuDisplay_1.OpacityArrayName = ['POINTS', 'Analytic']
s1_M5_PML_SHBC_1vtuDisplay_1.SelectInputVectors = ['POINTS', 'Displacement']
s1_M5_PML_SHBC_1vtuDisplay_1.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
s1_M5_PML_SHBC_1vtuDisplay_1.OSPRayScaleFunction.Points = [-0.0155398497892645, 0.0, 0.5, 0.0, -0.0009743210460688503, 0.5, 0.5, 0.0, 0.0135912076971268, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
s1_M5_PML_SHBC_1vtuDisplay_1.ScaleTransferFunction.Points = [-1.40464062906829, 0.0, 0.5, 0.0, 0.2805643188464151, 0.5, 0.5, 0.0, 1.9657692667611202, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
s1_M5_PML_SHBC_1vtuDisplay_1.OpacityTransferFunction.Points = [-1.40464062906829, 0.0, 0.5, 0.0, 0.2805643188464151, 0.5, 0.5, 0.0, 1.9657692667611202, 1.0, 0.5, 0.0]

# reset view to fit data
renderView2.ResetCamera()

# update the view to ensure updated data information
renderView2.Update()

# set active source
SetActiveSource(s1_M5_PML_SHBC_1_meshvtu)

# show data in view
s1_M5_PML_SHBC_1_meshvtuDisplay_1 = Show(s1_M5_PML_SHBC_1_meshvtu, renderView2, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
s1_M5_PML_SHBC_1_meshvtuDisplay_1.Representation = 'Surface'
s1_M5_PML_SHBC_1_meshvtuDisplay_1.ColorArrayName = [None, '']
s1_M5_PML_SHBC_1_meshvtuDisplay_1.SelectTCoordArray = 'None'
s1_M5_PML_SHBC_1_meshvtuDisplay_1.SelectNormalArray = 'None'
s1_M5_PML_SHBC_1_meshvtuDisplay_1.SelectTangentArray = 'None'
s1_M5_PML_SHBC_1_meshvtuDisplay_1.OSPRayScaleArray = 'Displacement'
s1_M5_PML_SHBC_1_meshvtuDisplay_1.OSPRayScaleFunction = 'PiecewiseFunction'
s1_M5_PML_SHBC_1_meshvtuDisplay_1.SelectOrientationVectors = 'Displacement'
s1_M5_PML_SHBC_1_meshvtuDisplay_1.ScaleFactor = 0.30000000000000004
s1_M5_PML_SHBC_1_meshvtuDisplay_1.SelectScaleArray = 'Displacement'
s1_M5_PML_SHBC_1_meshvtuDisplay_1.GlyphType = 'Arrow'
s1_M5_PML_SHBC_1_meshvtuDisplay_1.GlyphTableIndexArray = 'Displacement'
s1_M5_PML_SHBC_1_meshvtuDisplay_1.GaussianRadius = 0.015
s1_M5_PML_SHBC_1_meshvtuDisplay_1.SetScaleArray = ['POINTS', 'Displacement']
s1_M5_PML_SHBC_1_meshvtuDisplay_1.ScaleTransferFunction = 'PiecewiseFunction'
s1_M5_PML_SHBC_1_meshvtuDisplay_1.OpacityArray = ['POINTS', 'Displacement']
s1_M5_PML_SHBC_1_meshvtuDisplay_1.OpacityTransferFunction = 'PiecewiseFunction'
s1_M5_PML_SHBC_1_meshvtuDisplay_1.DataAxesGrid = 'GridAxesRepresentation'
s1_M5_PML_SHBC_1_meshvtuDisplay_1.PolarAxes = 'PolarAxesRepresentation'
s1_M5_PML_SHBC_1_meshvtuDisplay_1.ScalarOpacityUnitDistance = 0.2228352968694388
s1_M5_PML_SHBC_1_meshvtuDisplay_1.OpacityArrayName = ['POINTS', 'Displacement']
s1_M5_PML_SHBC_1_meshvtuDisplay_1.SelectInputVectors = ['POINTS', 'Displacement']
s1_M5_PML_SHBC_1_meshvtuDisplay_1.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
s1_M5_PML_SHBC_1_meshvtuDisplay_1.OSPRayScaleFunction.Points = [-0.0155398497892645, 0.0, 0.5, 0.0, -0.0009743210460688503, 0.5, 0.5, 0.0, 0.0135912076971268, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
s1_M5_PML_SHBC_1_meshvtuDisplay_1.ScaleTransferFunction.Points = [-4.08875818923397e-11, 0.0, 0.5, 0.0, 1.4187224806502317e-14, 0.5, 0.5, 0.0, 4.09159563419527e-11, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
s1_M5_PML_SHBC_1_meshvtuDisplay_1.OpacityTransferFunction.Points = [-4.08875818923397e-11, 0.0, 0.5, 0.0, 1.4187224806502317e-14, 0.5, 0.5, 0.0, 4.09159563419527e-11, 1.0, 0.5, 0.0]

# update the view to ensure updated data information
renderView2.Update()

# create a new 'Clip'
clip1 = Clip(registrationName='Clip1', Input=s1_M5_PML_SHBC_1_meshvtu)
clip1.ClipType = 'Plane'
clip1.HyperTreeGridClipper = 'Plane'
clip1.Scalars = [None, '']

# init the 'Plane' selected for 'ClipType'
clip1.ClipType.Origin = [0.0005201109570039941, 0.0, 0.0]

# init the 'Plane' selected for 'HyperTreeGridClipper'
clip1.HyperTreeGridClipper.Origin = [0.0005201109570039941, 0.0, 0.0]

# Properties modified on clip1
clip1.Scalars = ['POINTS', '']

# show data in view
clip1Display = Show(clip1, renderView2, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
clip1Display.Representation = 'Surface'
clip1Display.ColorArrayName = [None, '']
clip1Display.SelectTCoordArray = 'None'
clip1Display.SelectNormalArray = 'None'
clip1Display.SelectTangentArray = 'None'
clip1Display.OSPRayScaleArray = 'Displacement'
clip1Display.OSPRayScaleFunction = 'PiecewiseFunction'
clip1Display.SelectOrientationVectors = 'Displacement'
clip1Display.ScaleFactor = 0.30000000000000004
clip1Display.SelectScaleArray = 'Displacement'
clip1Display.GlyphType = 'Arrow'
clip1Display.GlyphTableIndexArray = 'Displacement'
clip1Display.GaussianRadius = 0.015
clip1Display.SetScaleArray = ['POINTS', 'Displacement']
clip1Display.ScaleTransferFunction = 'PiecewiseFunction'
clip1Display.OpacityArray = ['POINTS', 'Displacement']
clip1Display.OpacityTransferFunction = 'PiecewiseFunction'
clip1Display.DataAxesGrid = 'GridAxesRepresentation'
clip1Display.PolarAxes = 'PolarAxesRepresentation'
clip1Display.ScalarOpacityUnitDistance = 0.10609092398245626
clip1Display.OpacityArrayName = ['POINTS', 'Displacement']
clip1Display.SelectInputVectors = ['POINTS', 'Displacement']
clip1Display.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
clip1Display.OSPRayScaleFunction.Points = [-0.0155398497892645, 0.0, 0.5, 0.0, -0.0009743210460688503, 0.5, 0.5, 0.0, 0.0135912076971268, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
clip1Display.ScaleTransferFunction.Points = [-4.08875818923397e-11, 0.0, 0.5, 0.0, -5.263306788861799e-12, 0.5, 0.5, 0.0, 3.03609683146161e-11, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
clip1Display.OpacityTransferFunction.Points = [-4.08875818923397e-11, 0.0, 0.5, 0.0, -5.263306788861799e-12, 0.5, 0.5, 0.0, 3.03609683146161e-11, 1.0, 0.5, 0.0]

# hide data in view
Hide(s1_M5_PML_SHBC_1_meshvtu, renderView2)

# update the view to ensure updated data information
renderView2.Update()

# Properties modified on clip1.ClipType
clip1.ClipType.Normal = [1.0, 1.0, 0.0]

# update the view to ensure updated data information
renderView2.Update()

# Properties modified on clip1.ClipType
clip1.ClipType.Normal = [0.0, 1.0, 0.0]

# update the view to ensure updated data information
renderView2.Update()

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=clip1.ClipType)

# set active source
SetActiveSource(s1_M5_PML_SHBC_1vtu)

# create a new 'Clip'
clip2 = Clip(registrationName='Clip2', Input=s1_M5_PML_SHBC_1vtu)
clip2.ClipType = 'Plane'
clip2.HyperTreeGridClipper = 'Plane'
clip2.Scalars = ['POINTS', 'Analytic']
clip2.Value = 0.280564318846415

# init the 'Plane' selected for 'ClipType'
clip2.ClipType.Origin = [0.0005201109570039941, 0.0, 0.0]

# init the 'Plane' selected for 'HyperTreeGridClipper'
clip2.HyperTreeGridClipper.Origin = [0.0005201109570039941, 0.0, 0.0]

# show data in view
clip2Display = Show(clip2, renderView2, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
clip2Display.Representation = 'Surface'
clip2Display.ColorArrayName = [None, '']
clip2Display.SelectTCoordArray = 'None'
clip2Display.SelectNormalArray = 'None'
clip2Display.SelectTangentArray = 'None'
clip2Display.OSPRayScaleArray = 'Analytic'
clip2Display.OSPRayScaleFunction = 'PiecewiseFunction'
clip2Display.SelectOrientationVectors = 'Analytic'
clip2Display.ScaleFactor = 0.30000000000000004
clip2Display.SelectScaleArray = 'Analytic'
clip2Display.GlyphType = 'Arrow'
clip2Display.GlyphTableIndexArray = 'Analytic'
clip2Display.GaussianRadius = 0.015
clip2Display.SetScaleArray = ['POINTS', 'Analytic']
clip2Display.ScaleTransferFunction = 'PiecewiseFunction'
clip2Display.OpacityArray = ['POINTS', 'Analytic']
clip2Display.OpacityTransferFunction = 'PiecewiseFunction'
clip2Display.DataAxesGrid = 'GridAxesRepresentation'
clip2Display.PolarAxes = 'PolarAxesRepresentation'
clip2Display.ScalarOpacityUnitDistance = 0.13689161644163006
clip2Display.OpacityArrayName = ['POINTS', 'Analytic']
clip2Display.SelectInputVectors = ['POINTS', 'Displacement']
clip2Display.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
clip2Display.OSPRayScaleFunction.Points = [-0.0155398497892645, 0.0, 0.5, 0.0, -0.0009743210460688503, 0.5, 0.5, 0.0, 0.0135912076971268, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
clip2Display.ScaleTransferFunction.Points = [-1.40464062906829, 0.0, 0.5, 0.0, 0.2805643188464151, 0.5, 0.5, 0.0, 1.9657692667611202, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
clip2Display.OpacityTransferFunction.Points = [-1.40464062906829, 0.0, 0.5, 0.0, 0.2805643188464151, 0.5, 0.5, 0.0, 1.9657692667611202, 1.0, 0.5, 0.0]

# hide data in view
Hide(s1_M5_PML_SHBC_1vtu, renderView2)

# update the view to ensure updated data information
renderView2.Update()

# Properties modified on clip2.ClipType
clip2.ClipType.Normal = [1.0, 1.0, 0.0]

# update the view to ensure updated data information
renderView2.Update()

# Properties modified on clip2.ClipType
clip2.ClipType.Normal = [0.0, 1.0, 0.0]

# update the view to ensure updated data information
renderView2.Update()

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=clip2.ClipType)

# set active view
SetActiveView(renderView1)

# destroy renderView1
Delete(renderView1)
del renderView1

# close an empty frame
layout1.Collapse(1)

# set active view
SetActiveView(renderView2)

# create a new 'Clip'
clip3 = Clip(registrationName='Clip3', Input=clip2)
clip3.ClipType = 'Plane'
clip3.HyperTreeGridClipper = 'Plane'
clip3.Scalars = ['POINTS', 'Analytic']
clip3.Value = 0.280564318846415

# init the 'Plane' selected for 'ClipType'
clip3.ClipType.Origin = [0.0005201109570039941, -0.7498308756705035, 0.0]

# init the 'Plane' selected for 'HyperTreeGridClipper'
clip3.HyperTreeGridClipper.Origin = [0.0005201109570039941, -0.7498308756705035, 0.0]

# show data in view
clip3Display = Show(clip3, renderView2, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
clip3Display.Representation = 'Surface'
clip3Display.ColorArrayName = [None, '']
clip3Display.SelectTCoordArray = 'None'
clip3Display.SelectNormalArray = 'None'
clip3Display.SelectTangentArray = 'None'
clip3Display.OSPRayScaleArray = 'Analytic'
clip3Display.OSPRayScaleFunction = 'PiecewiseFunction'
clip3Display.SelectOrientationVectors = 'Analytic'
clip3Display.ScaleFactor = 0.30000000000000004
clip3Display.SelectScaleArray = 'Analytic'
clip3Display.GlyphType = 'Arrow'
clip3Display.GlyphTableIndexArray = 'Analytic'
clip3Display.GaussianRadius = 0.015
clip3Display.SetScaleArray = ['POINTS', 'Analytic']
clip3Display.ScaleTransferFunction = 'PiecewiseFunction'
clip3Display.OpacityArray = ['POINTS', 'Analytic']
clip3Display.OpacityTransferFunction = 'PiecewiseFunction'
clip3Display.DataAxesGrid = 'GridAxesRepresentation'
clip3Display.PolarAxes = 'PolarAxesRepresentation'
clip3Display.ScalarOpacityUnitDistance = 0.13935376454170872
clip3Display.OpacityArrayName = ['POINTS', 'Analytic']
clip3Display.SelectInputVectors = ['POINTS', 'Displacement']
clip3Display.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
clip3Display.OSPRayScaleFunction.Points = [-0.0155398497892645, 0.0, 0.5, 0.0, -0.0009743210460688503, 0.5, 0.5, 0.0, 0.0135912076971268, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
clip3Display.ScaleTransferFunction.Points = [-1.40464062906829, 0.0, 0.5, 0.0, 0.2805643188464151, 0.5, 0.5, 0.0, 1.9657692667611202, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
clip3Display.OpacityTransferFunction.Points = [-1.40464062906829, 0.0, 0.5, 0.0, 0.2805643188464151, 0.5, 0.5, 0.0, 1.9657692667611202, 1.0, 0.5, 0.0]

# hide data in view
Hide(clip2, renderView2)

# update the view to ensure updated data information
renderView2.Update()

# Properties modified on clip3.ClipType
clip3.ClipType.Normal = [-1.0, 0.0, 0.0]

# update the view to ensure updated data information
renderView2.Update()

# set active source
SetActiveSource(clip1)

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=clip3.ClipType)

# update the view to ensure updated data information
renderView2.Update()

# create a new 'Clip'
clip4 = Clip(registrationName='Clip4', Input=clip1)
clip4.ClipType = 'Plane'
clip4.HyperTreeGridClipper = 'Plane'
clip4.Scalars = [None, '']

# init the 'Plane' selected for 'ClipType'
clip4.ClipType.Origin = [0.0005201109570039941, -0.7498308756705023, 0.0]

# init the 'Plane' selected for 'HyperTreeGridClipper'
clip4.HyperTreeGridClipper.Origin = [0.0005201109570039941, -0.7498308756705023, 0.0]

# Properties modified on clip4
clip4.Scalars = ['POINTS', '']

# show data in view
clip4Display = Show(clip4, renderView2, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
clip4Display.Representation = 'Surface'
clip4Display.ColorArrayName = [None, '']
clip4Display.SelectTCoordArray = 'None'
clip4Display.SelectNormalArray = 'None'
clip4Display.SelectTangentArray = 'None'
clip4Display.OSPRayScaleArray = 'Displacement'
clip4Display.OSPRayScaleFunction = 'PiecewiseFunction'
clip4Display.SelectOrientationVectors = 'Displacement'
clip4Display.ScaleFactor = 0.30000000000000004
clip4Display.SelectScaleArray = 'Displacement'
clip4Display.GlyphType = 'Arrow'
clip4Display.GlyphTableIndexArray = 'Displacement'
clip4Display.GaussianRadius = 0.015
clip4Display.SetScaleArray = ['POINTS', 'Displacement']
clip4Display.ScaleTransferFunction = 'PiecewiseFunction'
clip4Display.OpacityArray = ['POINTS', 'Displacement']
clip4Display.OpacityTransferFunction = 'PiecewiseFunction'
clip4Display.DataAxesGrid = 'GridAxesRepresentation'
clip4Display.PolarAxes = 'PolarAxesRepresentation'
clip4Display.ScalarOpacityUnitDistance = 0.10817481473170923
clip4Display.OpacityArrayName = ['POINTS', 'Displacement']
clip4Display.SelectInputVectors = ['POINTS', 'Displacement']
clip4Display.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
clip4Display.OSPRayScaleFunction.Points = [-0.0155398497892645, 0.0, 0.5, 0.0, -0.0009743210460688503, 0.5, 0.5, 0.0, 0.0135912076971268, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
clip4Display.ScaleTransferFunction.Points = [-4.088758189233965e-11, 0.0, 0.5, 0.0, -5.2633067888617795e-12, 0.5, 0.5, 0.0, 3.036096831461609e-11, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
clip4Display.OpacityTransferFunction.Points = [-4.088758189233965e-11, 0.0, 0.5, 0.0, -5.2633067888617795e-12, 0.5, 0.5, 0.0, 3.036096831461609e-11, 1.0, 0.5, 0.0]

# hide data in view
Hide(clip1, renderView2)

# update the view to ensure updated data information
renderView2.Update()

# Properties modified on clip4.ClipType
clip4.ClipType.Normal = [-1.0, 0.0, 0.0]

# update the view to ensure updated data information
renderView2.Update()

# reset view to fit data
renderView2.ResetCamera()

# reset view to fit data
renderView2.ResetCamera()

#change interaction mode for render view
renderView2.InteractionMode = '2D'

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=clip4.ClipType)

# split cell
layout1.SplitVertical(0, 0.5)

# set active view
SetActiveView(None)

# set active view
SetActiveView(renderView2)

# set active source
SetActiveSource(clip3)

# toggle 3D widget visibility (only when running from the GUI)
Show3DWidgets(proxy=clip3.ClipType)

# update the view to ensure updated data information
renderView2.Update()

# set scalar coloring
ColorBy(clip3Display, ('POINTS', 'Scalar field'))

# rescale color and/or opacity maps used to include current data range
clip3Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
clip3Display.SetScalarBarVisibility(renderView2, True)

# get color transfer function/color map for 'Scalarfield'
scalarfieldLUT = GetColorTransferFunction('Scalarfield')

# get opacity transfer function/opacity map for 'Scalarfield'
scalarfieldPWF = GetOpacityTransferFunction('Scalarfield')

# set active view
SetActiveView(None)

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraFocalDisk = 1.0
# uncomment following to set a specific view size
# renderView1.ViewSize = [400, 400]

# assign view to a particular cell in the layout
AssignViewToLayout(view=renderView1, layout=layout1, hint=2)

# set active source
SetActiveSource(s1_M5_IE_SHBC_1vtu)

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=clip3.ClipType)

# update the view to ensure updated data information
renderView2.Update()

# create a new 'Clip'
clip5 = Clip(registrationName='Clip5', Input=s1_M5_IE_SHBC_1vtu)
clip5.ClipType = 'Plane'
clip5.HyperTreeGridClipper = 'Plane'
clip5.Scalars = ['POINTS', 'Analytic']
clip5.Value = 0.321408929335805

# init the 'Plane' selected for 'ClipType'
clip5.ClipType.Origin = [0.0005201109570039941, 0.0, 0.0]

# init the 'Plane' selected for 'HyperTreeGridClipper'
clip5.HyperTreeGridClipper.Origin = [0.0005201109570039941, 0.0, 0.0]

# show data in view
clip5Display = Show(clip5, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
clip5Display.Representation = 'Surface'
clip5Display.ColorArrayName = [None, '']
clip5Display.SelectTCoordArray = 'None'
clip5Display.SelectNormalArray = 'None'
clip5Display.SelectTangentArray = 'None'
clip5Display.OSPRayScaleArray = 'Analytic'
clip5Display.OSPRayScaleFunction = 'PiecewiseFunction'
clip5Display.SelectOrientationVectors = 'Analytic'
clip5Display.ScaleFactor = 0.30000000000000004
clip5Display.SelectScaleArray = 'Analytic'
clip5Display.GlyphType = 'Arrow'
clip5Display.GlyphTableIndexArray = 'Analytic'
clip5Display.GaussianRadius = 0.015
clip5Display.SetScaleArray = ['POINTS', 'Analytic']
clip5Display.ScaleTransferFunction = 'PiecewiseFunction'
clip5Display.OpacityArray = ['POINTS', 'Analytic']
clip5Display.OpacityTransferFunction = 'PiecewiseFunction'
clip5Display.DataAxesGrid = 'GridAxesRepresentation'
clip5Display.PolarAxes = 'PolarAxesRepresentation'
clip5Display.ScalarOpacityUnitDistance = 0.19178784665321952
clip5Display.OpacityArrayName = ['POINTS', 'Analytic']
clip5Display.SelectInputVectors = ['POINTS', 'Displacement']
clip5Display.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
clip5Display.OSPRayScaleFunction.Points = [-0.0155398497892645, 0.0, 0.5, 0.0, -0.0009743210460688503, 0.5, 0.5, 0.0, 0.0135912076971268, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
clip5Display.ScaleTransferFunction.Points = [-1.32295140808951, 0.0, 0.5, 0.0, 0.3214089293358051, 0.5, 0.5, 0.0, 1.9657692667611202, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
clip5Display.OpacityTransferFunction.Points = [-1.32295140808951, 0.0, 0.5, 0.0, 0.3214089293358051, 0.5, 0.5, 0.0, 1.9657692667611202, 1.0, 0.5, 0.0]

# reset view to fit data
renderView1.ResetCamera()

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on clip5.ClipType
clip5.ClipType.Normal = [1.0, 1.0, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on clip5.ClipType
clip5.ClipType.Normal = [0.0, 1.0, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Clip'
clip6 = Clip(registrationName='Clip6', Input=clip5)
clip6.ClipType = 'Plane'
clip6.HyperTreeGridClipper = 'Plane'
clip6.Scalars = ['POINTS', 'Analytic']
clip6.Value = 0.321408929335805

# init the 'Plane' selected for 'ClipType'
clip6.ClipType.Origin = [0.0005201109570039941, -0.7498308756705035, 0.0]

# init the 'Plane' selected for 'HyperTreeGridClipper'
clip6.HyperTreeGridClipper.Origin = [0.0005201109570039941, -0.7498308756705035, 0.0]

# show data in view
clip6Display = Show(clip6, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
clip6Display.Representation = 'Surface'
clip6Display.ColorArrayName = [None, '']
clip6Display.SelectTCoordArray = 'None'
clip6Display.SelectNormalArray = 'None'
clip6Display.SelectTangentArray = 'None'
clip6Display.OSPRayScaleArray = 'Analytic'
clip6Display.OSPRayScaleFunction = 'PiecewiseFunction'
clip6Display.SelectOrientationVectors = 'Analytic'
clip6Display.ScaleFactor = 0.30000000000000004
clip6Display.SelectScaleArray = 'Analytic'
clip6Display.GlyphType = 'Arrow'
clip6Display.GlyphTableIndexArray = 'Analytic'
clip6Display.GaussianRadius = 0.015
clip6Display.SetScaleArray = ['POINTS', 'Analytic']
clip6Display.ScaleTransferFunction = 'PiecewiseFunction'
clip6Display.OpacityArray = ['POINTS', 'Analytic']
clip6Display.OpacityTransferFunction = 'PiecewiseFunction'
clip6Display.DataAxesGrid = 'GridAxesRepresentation'
clip6Display.PolarAxes = 'PolarAxesRepresentation'
clip6Display.ScalarOpacityUnitDistance = 0.19442639768466188
clip6Display.OpacityArrayName = ['POINTS', 'Analytic']
clip6Display.SelectInputVectors = ['POINTS', 'Displacement']
clip6Display.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
clip6Display.OSPRayScaleFunction.Points = [-0.0155398497892645, 0.0, 0.5, 0.0, -0.0009743210460688503, 0.5, 0.5, 0.0, 0.0135912076971268, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
clip6Display.ScaleTransferFunction.Points = [-1.32295140808951, 0.0, 0.5, 0.0, 0.3214089293358051, 0.5, 0.5, 0.0, 1.9657692667611202, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
clip6Display.OpacityTransferFunction.Points = [-1.32295140808951, 0.0, 0.5, 0.0, 0.3214089293358051, 0.5, 0.5, 0.0, 1.9657692667611202, 1.0, 0.5, 0.0]

# hide data in view
Hide(clip5, renderView1)

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on clip6.ClipType
clip6.ClipType.Normal = [-1.0, 0.0, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# reset view to fit data
renderView1.ResetCamera()

# reset view to fit data
renderView1.ResetCamera()

# set active view
SetActiveView(renderView2)

#change interaction mode for render view
renderView2.InteractionMode = '3D'

#change interaction mode for render view
renderView2.InteractionMode = '2D'

# set active view
SetActiveView(renderView1)

#change interaction mode for render view
renderView1.InteractionMode = '2D'

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=clip6.ClipType)

# set active source
SetActiveSource(s1_M5_IE_SHBC_1_meshvtu)

# update the view to ensure updated data information
renderView2.Update()

# create a new 'Clip'
clip7 = Clip(registrationName='Clip7', Input=s1_M5_IE_SHBC_1_meshvtu)
clip7.ClipType = 'Plane'
clip7.HyperTreeGridClipper = 'Plane'
clip7.Scalars = [None, '']

# init the 'Plane' selected for 'ClipType'
clip7.ClipType.Origin = [0.0005201109570039941, 0.0, 0.0]

# init the 'Plane' selected for 'HyperTreeGridClipper'
clip7.HyperTreeGridClipper.Origin = [0.0005201109570039941, 0.0, 0.0]

# Properties modified on clip7
clip7.Scalars = ['POINTS', '']

# show data in view
clip7Display = Show(clip7, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
clip7Display.Representation = 'Surface'
clip7Display.ColorArrayName = [None, '']
clip7Display.SelectTCoordArray = 'None'
clip7Display.SelectNormalArray = 'None'
clip7Display.SelectTangentArray = 'None'
clip7Display.OSPRayScaleArray = 'Displacement'
clip7Display.OSPRayScaleFunction = 'PiecewiseFunction'
clip7Display.SelectOrientationVectors = 'Displacement'
clip7Display.ScaleFactor = 0.30000000000000004
clip7Display.SelectScaleArray = 'Displacement'
clip7Display.GlyphType = 'Arrow'
clip7Display.GlyphTableIndexArray = 'Displacement'
clip7Display.GaussianRadius = 0.015
clip7Display.SetScaleArray = ['POINTS', 'Displacement']
clip7Display.ScaleTransferFunction = 'PiecewiseFunction'
clip7Display.OpacityArray = ['POINTS', 'Displacement']
clip7Display.OpacityTransferFunction = 'PiecewiseFunction'
clip7Display.DataAxesGrid = 'GridAxesRepresentation'
clip7Display.PolarAxes = 'PolarAxesRepresentation'
clip7Display.ScalarOpacityUnitDistance = 0.1490971963936087
clip7Display.OpacityArrayName = ['POINTS', 'Displacement']
clip7Display.SelectInputVectors = ['POINTS', 'Displacement']
clip7Display.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
clip7Display.OSPRayScaleFunction.Points = [-0.0155398497892645, 0.0, 0.5, 0.0, -0.0009743210460688503, 0.5, 0.5, 0.0, 0.0135912076971268, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
clip7Display.ScaleTransferFunction.Points = [-4.13450632764285e-11, 0.0, 0.5, 0.0, -4.7013250744076e-12, 0.5, 0.5, 0.0, 3.19424131276133e-11, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
clip7Display.OpacityTransferFunction.Points = [-4.13450632764285e-11, 0.0, 0.5, 0.0, -4.7013250744076e-12, 0.5, 0.5, 0.0, 3.19424131276133e-11, 1.0, 0.5, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on clip7.ClipType
clip7.ClipType.Normal = [1.0, 1.0, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Clip'
clip8 = Clip(registrationName='Clip8', Input=clip7)
clip8.ClipType = 'Plane'
clip8.HyperTreeGridClipper = 'Plane'
clip8.Scalars = [None, '']

# init the 'Plane' selected for 'ClipType'
clip8.ClipType.Origin = [-0.21907126546703526, -0.21962301009363194, 0.0]

# init the 'Plane' selected for 'HyperTreeGridClipper'
clip8.HyperTreeGridClipper.Origin = [-0.21907126546703526, -0.21962301009363194, 0.0]

# Properties modified on clip7.ClipType
clip7.ClipType.Normal = [0.0, 1.0, 0.0]

# Properties modified on clip8
clip8.Scalars = ['POINTS', '']

# show data in view
clip8Display = Show(clip8, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
clip8Display.Representation = 'Surface'
clip8Display.ColorArrayName = [None, '']
clip8Display.SelectTCoordArray = 'None'
clip8Display.SelectNormalArray = 'None'
clip8Display.SelectTangentArray = 'None'
clip8Display.OSPRayScaleArray = 'Displacement'
clip8Display.OSPRayScaleFunction = 'PiecewiseFunction'
clip8Display.SelectOrientationVectors = 'Displacement'
clip8Display.ScaleFactor = 0.2949916105249724
clip8Display.SelectScaleArray = 'Displacement'
clip8Display.GlyphType = 'Arrow'
clip8Display.GlyphTableIndexArray = 'Displacement'
clip8Display.GaussianRadius = 0.014749580526248621
clip8Display.SetScaleArray = ['POINTS', 'Displacement']
clip8Display.ScaleTransferFunction = 'PiecewiseFunction'
clip8Display.OpacityArray = ['POINTS', 'Displacement']
clip8Display.OpacityTransferFunction = 'PiecewiseFunction'
clip8Display.DataAxesGrid = 'GridAxesRepresentation'
clip8Display.PolarAxes = 'PolarAxesRepresentation'
clip8Display.ScalarOpacityUnitDistance = 0.163817746050322
clip8Display.OpacityArrayName = ['POINTS', 'Displacement']
clip8Display.SelectInputVectors = ['POINTS', 'Displacement']
clip8Display.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
clip8Display.OSPRayScaleFunction.Points = [-0.0155398497892645, 0.0, 0.5, 0.0, -0.0009743210460688503, 0.5, 0.5, 0.0, 0.0135912076971268, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
clip8Display.ScaleTransferFunction.Points = [-4.13450632764285e-11, 0.0, 0.5, 0.0, -4.7013250744076e-12, 0.5, 0.5, 0.0, 3.19424131276133e-11, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
clip8Display.OpacityTransferFunction.Points = [-4.13450632764285e-11, 0.0, 0.5, 0.0, -4.7013250744076e-12, 0.5, 0.5, 0.0, 3.19424131276133e-11, 1.0, 0.5, 0.0]

# hide data in view
Hide(clip7, renderView1)

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on clip8.ClipType
clip8.ClipType.Normal = [-1.0, 0.0, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(clip7)

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=clip8.ClipType)

# toggle 3D widget visibility (only when running from the GUI)
Show3DWidgets(proxy=clip7.ClipType)

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(clip8)

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=clip7.ClipType)

# toggle 3D widget visibility (only when running from the GUI)
Show3DWidgets(proxy=clip8.ClipType)

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(clip7)

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=clip8.ClipType)

# toggle 3D widget visibility (only when running from the GUI)
Show3DWidgets(proxy=clip7.ClipType)

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(clip8)

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=clip7.ClipType)

# toggle 3D widget visibility (only when running from the GUI)
Show3DWidgets(proxy=clip8.ClipType)

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on clip8.ClipType
clip8.ClipType.Origin = [0.0, -0.21962301009363194, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on clip8.ClipType
clip8.ClipType.Origin = [0.0, 0.0, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(clip7)

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=clip8.ClipType)

# toggle 3D widget visibility (only when running from the GUI)
Show3DWidgets(proxy=clip7.ClipType)

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on clip7.ClipType
clip7.ClipType.Origin = [0.0, 0.0, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(clip5)

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=clip7.ClipType)

# toggle 3D widget visibility (only when running from the GUI)
Show3DWidgets(proxy=clip5.ClipType)

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on clip5.ClipType
clip5.ClipType.Origin = [0.0, 0.0, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(s1_M5_IE_SHBC_1vtu)

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=clip5.ClipType)

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(clip5)

# toggle 3D widget visibility (only when running from the GUI)
Show3DWidgets(proxy=clip5.ClipType)

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(clip6)

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=clip5.ClipType)

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on clip6.ClipType
clip6.ClipType.Origin = [0.0, -0.7498308756705035, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on clip6.ClipType
clip6.ClipType.Origin = [0.0, 0.0, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(clip8)

# toggle 3D widget visibility (only when running from the GUI)
Show3DWidgets(proxy=clip8.ClipType)

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(clip6)

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=clip8.ClipType)

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(clip6Display, ('POINTS', 'Scalar field'))

# rescale color and/or opacity maps used to include current data range
clip6Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
clip6Display.SetScalarBarVisibility(renderView1, True)

# create a new 'Sphere'
sphere1 = Sphere(registrationName='Sphere1')

# show data in view
sphere1Display = Show(sphere1, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
sphere1Display.Representation = 'Surface'
sphere1Display.ColorArrayName = [None, '']
sphere1Display.SelectTCoordArray = 'None'
sphere1Display.SelectNormalArray = 'Normals'
sphere1Display.SelectTangentArray = 'None'
sphere1Display.OSPRayScaleArray = 'Normals'
sphere1Display.OSPRayScaleFunction = 'PiecewiseFunction'
sphere1Display.SelectOrientationVectors = 'None'
sphere1Display.ScaleFactor = 0.1
sphere1Display.SelectScaleArray = 'None'
sphere1Display.GlyphType = 'Arrow'
sphere1Display.GlyphTableIndexArray = 'None'
sphere1Display.GaussianRadius = 0.005
sphere1Display.SetScaleArray = ['POINTS', 'Normals']
sphere1Display.ScaleTransferFunction = 'PiecewiseFunction'
sphere1Display.OpacityArray = ['POINTS', 'Normals']
sphere1Display.OpacityTransferFunction = 'PiecewiseFunction'
sphere1Display.DataAxesGrid = 'GridAxesRepresentation'
sphere1Display.PolarAxes = 'PolarAxesRepresentation'
sphere1Display.SelectInputVectors = ['POINTS', 'Normals']
sphere1Display.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
sphere1Display.OSPRayScaleFunction.Points = [-0.0155398497892645, 0.0, 0.5, 0.0, -0.0009743210460688503, 0.5, 0.5, 0.0, 0.0135912076971268, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
sphere1Display.ScaleTransferFunction.Points = [-0.9749279022216797, 0.0, 0.5, 0.0, 0.0, 0.5, 0.5, 0.0, 0.9749279022216797, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
sphere1Display.OpacityTransferFunction.Points = [-0.9749279022216797, 0.0, 0.5, 0.0, 0.0, 0.5, 0.5, 0.0, 0.9749279022216797, 1.0, 0.5, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on sphere1
sphere1.Radius = 1.0

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on sphere1
sphere1.ThetaResolution = 100

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on sphere1
sphere1.PhiResolution = 100

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on sphere1
sphere1.ThetaResolution = 10

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on sphere1
sphere1.ThetaResolution = 2

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on sphere1
sphere1.ThetaResolution = 5

# update the view to ensure updated data information
renderView1.Update()

# change representation type
sphere1Display.SetRepresentationType('Surface With Edges')

# set active view
SetActiveView(renderView2)

#change interaction mode for render view
renderView2.InteractionMode = '3D'

#change interaction mode for render view
renderView2.InteractionMode = '2D'

# set active view
SetActiveView(renderView1)

#change interaction mode for render view
renderView1.InteractionMode = '3D'

#change interaction mode for render view
renderView1.InteractionMode = '2D'

#change interaction mode for render view
renderView1.InteractionMode = '3D'

# create a new 'Clip'
clip9 = Clip(registrationName='Clip9', Input=sphere1)
clip9.ClipType = 'Plane'
clip9.HyperTreeGridClipper = 'Plane'
clip9.Scalars = [None, '']

# init the 'Plane' selected for 'ClipType'
clip9.ClipType.Origin = [0.09547948837280273, 0.0, 0.0]

# init the 'Plane' selected for 'HyperTreeGridClipper'
clip9.HyperTreeGridClipper.Origin = [0.09547948837280273, 0.0, 0.0]

# Properties modified on clip9
clip9.Scalars = ['POINTS', '']

# show data in view
clip9Display = Show(clip9, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
clip9Display.Representation = 'Surface'
clip9Display.ColorArrayName = [None, '']
clip9Display.SelectTCoordArray = 'None'
clip9Display.SelectNormalArray = 'Normals'
clip9Display.SelectTangentArray = 'None'
clip9Display.OSPRayScaleArray = 'Normals'
clip9Display.OSPRayScaleFunction = 'PiecewiseFunction'
clip9Display.SelectOrientationVectors = 'None'
clip9Display.ScaleFactor = 0.2
clip9Display.SelectScaleArray = 'None'
clip9Display.GlyphType = 'Arrow'
clip9Display.GlyphTableIndexArray = 'None'
clip9Display.GaussianRadius = 0.01
clip9Display.SetScaleArray = ['POINTS', 'Normals']
clip9Display.ScaleTransferFunction = 'PiecewiseFunction'
clip9Display.OpacityArray = ['POINTS', 'Normals']
clip9Display.OpacityTransferFunction = 'PiecewiseFunction'
clip9Display.DataAxesGrid = 'GridAxesRepresentation'
clip9Display.PolarAxes = 'PolarAxesRepresentation'
clip9Display.ScalarOpacityUnitDistance = 0.3230402678360829
clip9Display.OpacityArrayName = ['POINTS', 'Normals']
clip9Display.SelectInputVectors = ['POINTS', 'Normals']
clip9Display.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
clip9Display.OSPRayScaleFunction.Points = [-0.0155398497892645, 0.0, 0.5, 0.0, -0.0009743210460688503, 0.5, 0.5, 0.0, 0.0135912076971268, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
clip9Display.ScaleTransferFunction.Points = [-0.8089151382446289, 0.0, 0.5, 0.0, -0.3567178249359131, 0.5, 0.5, 0.0, 0.09547948837280273, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
clip9Display.OpacityTransferFunction.Points = [-0.8089151382446289, 0.0, 0.5, 0.0, -0.3567178249359131, 0.5, 0.5, 0.0, 0.09547948837280273, 1.0, 0.5, 0.0]

# hide data in view
Hide(sphere1, renderView1)

# update the view to ensure updated data information
renderView2.Update()

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on clip9.ClipType
clip9.ClipType.Normal = [-1.0, 0.0, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Slice'
slice1 = Slice(registrationName='Slice1', Input=clip9)
slice1.SliceType = 'Plane'
slice1.HyperTreeGridSlicer = 'Plane'
slice1.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1.SliceType.Origin = [0.5476768016815186, 0.0, 0.0]

# init the 'Plane' selected for 'HyperTreeGridSlicer'
slice1.HyperTreeGridSlicer.Origin = [0.5476768016815186, 0.0, 0.0]

# show data in view
slice1Display = Show(slice1, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
slice1Display.Representation = 'Surface'
slice1Display.ColorArrayName = [None, '']
slice1Display.SelectTCoordArray = 'None'
slice1Display.SelectNormalArray = 'Normals'
slice1Display.SelectTangentArray = 'None'
slice1Display.OSPRayScaleArray = 'Normals'
slice1Display.OSPRayScaleFunction = 'PiecewiseFunction'
slice1Display.SelectOrientationVectors = 'None'
slice1Display.ScaleFactor = 0.1673144817352295
slice1Display.SelectScaleArray = 'None'
slice1Display.GlyphType = 'Arrow'
slice1Display.GlyphTableIndexArray = 'None'
slice1Display.GaussianRadius = 0.008365724086761474
slice1Display.SetScaleArray = ['POINTS', 'Normals']
slice1Display.ScaleTransferFunction = 'PiecewiseFunction'
slice1Display.OpacityArray = ['POINTS', 'Normals']
slice1Display.OpacityTransferFunction = 'PiecewiseFunction'
slice1Display.DataAxesGrid = 'GridAxesRepresentation'
slice1Display.PolarAxes = 'PolarAxesRepresentation'
slice1Display.SelectInputVectors = ['POINTS', 'Normals']
slice1Display.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
slice1Display.OSPRayScaleFunction.Points = [-0.0155398497892645, 0.0, 0.5, 0.0, -0.0009743210460688503, 0.5, 0.5, 0.0, 0.0135912076971268, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
slice1Display.ScaleTransferFunction.Points = [0.5476768016815186, 0.0, 0.5, 0.0, 0.5477378368377686, 0.5, 0.5, 0.0, 0.5477988719940186, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
slice1Display.OpacityTransferFunction.Points = [0.5476768016815186, 0.0, 0.5, 0.0, 0.5477378368377686, 0.5, 0.5, 0.0, 0.5477988719940186, 1.0, 0.5, 0.0]

# hide data in view
Hide(clip9, renderView1)

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on slice1.SliceType
slice1.SliceType.Normal = [1.0, 1.0, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on slice1.SliceType
slice1.SliceType.Normal = [0.0, 1.0, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Resample With Dataset'
resampleWithDataset1 = ResampleWithDataset(registrationName='ResampleWithDataset1', SourceDataArrays=clip6,
    DestinationMesh=slice1)
resampleWithDataset1.CellLocator = 'Static Cell Locator'

# show data in view
resampleWithDataset1Display = Show(resampleWithDataset1, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
resampleWithDataset1Display.Representation = 'Surface'
resampleWithDataset1Display.ColorArrayName = [None, '']
resampleWithDataset1Display.SelectTCoordArray = 'None'
resampleWithDataset1Display.SelectNormalArray = 'None'
resampleWithDataset1Display.SelectTangentArray = 'None'
resampleWithDataset1Display.OSPRayScaleArray = 'Analytic'
resampleWithDataset1Display.OSPRayScaleFunction = 'PiecewiseFunction'
resampleWithDataset1Display.SelectOrientationVectors = 'Analytic'
resampleWithDataset1Display.ScaleFactor = 0.19908493757247925
resampleWithDataset1Display.SelectScaleArray = 'Analytic'
resampleWithDataset1Display.GlyphType = 'Arrow'
resampleWithDataset1Display.GlyphTableIndexArray = 'Analytic'
resampleWithDataset1Display.GaussianRadius = 0.009954246878623963
resampleWithDataset1Display.SetScaleArray = ['POINTS', 'Analytic']
resampleWithDataset1Display.ScaleTransferFunction = 'PiecewiseFunction'
resampleWithDataset1Display.OpacityArray = ['POINTS', 'Analytic']
resampleWithDataset1Display.OpacityTransferFunction = 'PiecewiseFunction'
resampleWithDataset1Display.DataAxesGrid = 'GridAxesRepresentation'
resampleWithDataset1Display.PolarAxes = 'PolarAxesRepresentation'
resampleWithDataset1Display.SelectInputVectors = ['POINTS', 'Displacement']
resampleWithDataset1Display.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
resampleWithDataset1Display.OSPRayScaleFunction.Points = [-0.0155398497892645, 0.0, 0.5, 0.0, -0.0009743210460688503, 0.5, 0.5, 0.0, 0.0135912076971268, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
resampleWithDataset1Display.ScaleTransferFunction.Points = [-0.9503156600884082, 0.0, 0.5, 0.0, 0.2530971762314853, 0.5, 0.5, 0.0, 1.4565100125513788, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
resampleWithDataset1Display.OpacityTransferFunction.Points = [-0.9503156600884082, 0.0, 0.5, 0.0, 0.2530971762314853, 0.5, 0.5, 0.0, 1.4565100125513788, 1.0, 0.5, 0.0]

# hide data in view
Hide(slice1, renderView1)

# hide data in view
Hide(clip6, renderView1)

# update the view to ensure updated data information
renderView1.Update()

CreateLayout('Layout #2')

# set active view
SetActiveView(None)

# create a new 'Plot Data'
plotData1 = PlotData(registrationName='PlotData1', Input=resampleWithDataset1)

# Create a new 'Line Chart View'
lineChartView1 = CreateView('XYChartView')
# uncomment following to set a specific view size
# lineChartView1.ViewSize = [400, 400]

# show data in view
plotData1Display = Show(plotData1, lineChartView1, 'XYChartRepresentation')

# trace defaults for the display properties.
plotData1Display.CompositeDataSetIndex = [0]
plotData1Display.XArrayName = 'Analytic'
plotData1Display.SeriesVisibility = ['Analytic', 'Displacement_Magnitude', 'Error', 'Error Gradient', 'Error in Energy', 'P_inc', 'Scalar field', 'Total scalar field (abs)', 'Total scalar field (real)']
plotData1Display.SeriesLabel = ['Analytic', 'Analytic', 'Displacement_X', 'Displacement_X', 'Displacement_Y', 'Displacement_Y', 'Displacement_Z', 'Displacement_Z', 'Displacement_Magnitude', 'Displacement_Magnitude', 'Error', 'Error', 'Error Gradient', 'Error Gradient', 'Error in Energy', 'Error in Energy', 'P_inc', 'P_inc', 'Scalar field', 'Scalar field', 'Total scalar field (abs)', 'Total scalar field (abs)', 'Total scalar field (real)', 'Total scalar field (real)', 'vtkValidPointMask', 'vtkValidPointMask', 'Points_X', 'Points_X', 'Points_Y', 'Points_Y', 'Points_Z', 'Points_Z', 'Points_Magnitude', 'Points_Magnitude']
plotData1Display.SeriesColor = ['Analytic', '0', '0', '0', 'Displacement_X', '0.8899977111467154', '0.10000762951094835', '0.1100022888532845', 'Displacement_Y', '0.220004577706569', '0.4899977111467155', '0.7199969481956207', 'Displacement_Z', '0.30000762951094834', '0.6899977111467155', '0.2899977111467155', 'Displacement_Magnitude', '0.6', '0.3100022888532845', '0.6399938963912413', 'Error', '1', '0.5000076295109483', '0', 'Error Gradient', '0.6500038147554742', '0.3400015259021897', '0.16000610360875867', 'Error in Energy', '0', '0', '0', 'P_inc', '0.8899977111467154', '0.10000762951094835', '0.1100022888532845', 'Scalar field', '0.220004577706569', '0.4899977111467155', '0.7199969481956207', 'Total scalar field (abs)', '0.30000762951094834', '0.6899977111467155', '0.2899977111467155', 'Total scalar field (real)', '0.6', '0.3100022888532845', '0.6399938963912413', 'vtkValidPointMask', '1', '0.5000076295109483', '0', 'Points_X', '0.6500038147554742', '0.3400015259021897', '0.16000610360875867', 'Points_Y', '0', '0', '0', 'Points_Z', '0.8899977111467154', '0.10000762951094835', '0.1100022888532845', 'Points_Magnitude', '0.220004577706569', '0.4899977111467155', '0.7199969481956207']
plotData1Display.SeriesPlotCorner = ['Analytic', '0', 'Displacement_X', '0', 'Displacement_Y', '0', 'Displacement_Z', '0', 'Displacement_Magnitude', '0', 'Error', '0', 'Error Gradient', '0', 'Error in Energy', '0', 'P_inc', '0', 'Scalar field', '0', 'Total scalar field (abs)', '0', 'Total scalar field (real)', '0', 'vtkValidPointMask', '0', 'Points_X', '0', 'Points_Y', '0', 'Points_Z', '0', 'Points_Magnitude', '0']
plotData1Display.SeriesLabelPrefix = ''
plotData1Display.SeriesLineStyle = ['Analytic', '1', 'Displacement_X', '1', 'Displacement_Y', '1', 'Displacement_Z', '1', 'Displacement_Magnitude', '1', 'Error', '1', 'Error Gradient', '1', 'Error in Energy', '1', 'P_inc', '1', 'Scalar field', '1', 'Total scalar field (abs)', '1', 'Total scalar field (real)', '1', 'vtkValidPointMask', '1', 'Points_X', '1', 'Points_Y', '1', 'Points_Z', '1', 'Points_Magnitude', '1']
plotData1Display.SeriesLineThickness = ['Analytic', '2', 'Displacement_X', '2', 'Displacement_Y', '2', 'Displacement_Z', '2', 'Displacement_Magnitude', '2', 'Error', '2', 'Error Gradient', '2', 'Error in Energy', '2', 'P_inc', '2', 'Scalar field', '2', 'Total scalar field (abs)', '2', 'Total scalar field (real)', '2', 'vtkValidPointMask', '2', 'Points_X', '2', 'Points_Y', '2', 'Points_Z', '2', 'Points_Magnitude', '2']
plotData1Display.SeriesMarkerStyle = ['Analytic', '0', 'Displacement_X', '0', 'Displacement_Y', '0', 'Displacement_Z', '0', 'Displacement_Magnitude', '0', 'Error', '0', 'Error Gradient', '0', 'Error in Energy', '0', 'P_inc', '0', 'Scalar field', '0', 'Total scalar field (abs)', '0', 'Total scalar field (real)', '0', 'vtkValidPointMask', '0', 'Points_X', '0', 'Points_Y', '0', 'Points_Z', '0', 'Points_Magnitude', '0']
plotData1Display.SeriesMarkerSize = ['Analytic', '4', 'Displacement_X', '4', 'Displacement_Y', '4', 'Displacement_Z', '4', 'Displacement_Magnitude', '4', 'Error', '4', 'Error Gradient', '4', 'Error in Energy', '4', 'P_inc', '4', 'Scalar field', '4', 'Total scalar field (abs)', '4', 'Total scalar field (real)', '4', 'vtkValidPointMask', '4', 'Points_X', '4', 'Points_Y', '4', 'Points_Z', '4', 'Points_Magnitude', '4']

# get layout
layout2 = GetLayoutByName("Layout #2")

# add view to a layout so it's visible in UI
AssignViewToLayout(view=lineChartView1, layout=layout2, hint=0)

# Properties modified on plotData1Display
plotData1Display.SeriesPlotCorner = ['Analytic', '0', 'Displacement_Magnitude', '0', 'Displacement_X', '0', 'Displacement_Y', '0', 'Displacement_Z', '0', 'Error', '0', 'Error Gradient', '0', 'Error in Energy', '0', 'P_inc', '0', 'Points_Magnitude', '0', 'Points_X', '0', 'Points_Y', '0', 'Points_Z', '0', 'Scalar field', '0', 'Total scalar field (abs)', '0', 'Total scalar field (real)', '0', 'vtkValidPointMask', '0']
plotData1Display.SeriesLineStyle = ['Analytic', '1', 'Displacement_Magnitude', '1', 'Displacement_X', '1', 'Displacement_Y', '1', 'Displacement_Z', '1', 'Error', '1', 'Error Gradient', '1', 'Error in Energy', '1', 'P_inc', '1', 'Points_Magnitude', '1', 'Points_X', '1', 'Points_Y', '1', 'Points_Z', '1', 'Scalar field', '1', 'Total scalar field (abs)', '1', 'Total scalar field (real)', '1', 'vtkValidPointMask', '1']
plotData1Display.SeriesLineThickness = ['Analytic', '2', 'Displacement_Magnitude', '2', 'Displacement_X', '2', 'Displacement_Y', '2', 'Displacement_Z', '2', 'Error', '2', 'Error Gradient', '2', 'Error in Energy', '2', 'P_inc', '2', 'Points_Magnitude', '2', 'Points_X', '2', 'Points_Y', '2', 'Points_Z', '2', 'Scalar field', '2', 'Total scalar field (abs)', '2', 'Total scalar field (real)', '2', 'vtkValidPointMask', '2']
plotData1Display.SeriesMarkerStyle = ['Analytic', '0', 'Displacement_Magnitude', '0', 'Displacement_X', '0', 'Displacement_Y', '0', 'Displacement_Z', '0', 'Error', '0', 'Error Gradient', '0', 'Error in Energy', '0', 'P_inc', '0', 'Points_Magnitude', '0', 'Points_X', '0', 'Points_Y', '0', 'Points_Z', '0', 'Scalar field', '0', 'Total scalar field (abs)', '0', 'Total scalar field (real)', '0', 'vtkValidPointMask', '0']
plotData1Display.SeriesMarkerSize = ['Analytic', '4', 'Displacement_Magnitude', '4', 'Displacement_X', '4', 'Displacement_Y', '4', 'Displacement_Z', '4', 'Error', '4', 'Error Gradient', '4', 'Error in Energy', '4', 'P_inc', '4', 'Points_Magnitude', '4', 'Points_X', '4', 'Points_Y', '4', 'Points_Z', '4', 'Scalar field', '4', 'Total scalar field (abs)', '4', 'Total scalar field (real)', '4', 'vtkValidPointMask', '4']

# update the view to ensure updated data information
lineChartView1.Update()

# set active view
SetActiveView(renderView2)

# set active source
SetActiveSource(resampleWithDataset1)

# set active view
SetActiveView(renderView1)

# set active source
SetActiveSource(clip6)

# show data in view
clip6Display = Show(clip6, renderView1, 'UnstructuredGridRepresentation')

# show color bar/color legend
clip6Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# hide data in view
Hide(clip6, renderView1)

# hide data in view
Hide(clip8, renderView1)

# set active source
SetActiveSource(clip8)

# toggle 3D widget visibility (only when running from the GUI)
Show3DWidgets(proxy=clip8.ClipType)

# show data in view
clip8Display = Show(clip8, renderView1, 'UnstructuredGridRepresentation')

# update the view to ensure updated data information
renderView1.Update()

# hide data in view
Hide(clip8, renderView1)

# set active source
SetActiveSource(clip6)

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=clip8.ClipType)

# show data in view
clip6Display = Show(clip6, renderView1, 'UnstructuredGridRepresentation')

# show color bar/color legend
clip6Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(sphere1)

# Properties modified on sphere1
sphere1.ThetaResolution = 100

# update the view to ensure updated data information
renderView1.Update()

# update the view to ensure updated data information
lineChartView1.Update()

# Properties modified on sphere1
sphere1.ThetaResolution = 5

# update the view to ensure updated data information
renderView1.Update()

# update the view to ensure updated data information
lineChartView1.Update()

# Properties modified on sphere1
sphere1.PhiResolution = 1000

# update the view to ensure updated data information
renderView1.Update()

# update the view to ensure updated data information
lineChartView1.Update()

# set active source
SetActiveSource(clip4)

# set active source
SetActiveSource(clip9)

# toggle 3D widget visibility (only when running from the GUI)
Show3DWidgets(proxy=clip9.ClipType)

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on clip9.ClipType
clip9.ClipType.Origin = [0.0, 0.0, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# update the view to ensure updated data information
lineChartView1.Update()

# set active source
SetActiveSource(slice1)

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=clip9.ClipType)

# toggle 3D widget visibility (only when running from the GUI)
Show3DWidgets(proxy=slice1.SliceType)

# Properties modified on slice1.SliceType
slice1.SliceType.Origin = [0.0, 0.0, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# update the view to ensure updated data information
lineChartView1.Update()

# set active source
SetActiveSource(resampleWithDataset1)

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=slice1.SliceType)

# set active view
SetActiveView(lineChartView1)

# set active source
SetActiveSource(plotData1)

# Properties modified on plotData1Display
plotData1Display.UseIndexForXAxis = 0

# set active view
SetActiveView(renderView2)

# set active view
SetActiveView(lineChartView1)

# Properties modified on plotData1Display
plotData1Display.XArrayName = 'Points_Z'

# set active source
SetActiveSource(resampleWithDataset1)

# update the view to ensure updated data information
renderView2.Update()

# set active view
SetActiveView(renderView2)

# set active source
SetActiveSource(plotData1)

# set active source
SetActiveSource(resampleWithDataset1)

# set active view
SetActiveView(renderView1)

# set scalar coloring
ColorBy(resampleWithDataset1Display, ('POINTS', 'Scalar field'))

# rescale color and/or opacity maps used to include current data range
resampleWithDataset1Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
resampleWithDataset1Display.SetScalarBarVisibility(renderView1, True)

# set active source
SetActiveSource(sphere1)

# Properties modified on sphere1
sphere1.Radius = 1.01

# update the view to ensure updated data information
renderView1.Update()

# update the view to ensure updated data information
lineChartView1.Update()

# Properties modified on sphere1
sphere1.ThetaResolution = 100

# update the view to ensure updated data information
renderView1.Update()

# update the view to ensure updated data information
lineChartView1.Update()

# set active view
SetActiveView(lineChartView1)

# set active view
SetActiveView(renderView2)

# set active view
SetActiveView(renderView1)

# hide data in view
Hide(clip6, renderView1)

# set active source
SetActiveSource(clip6)

# show data in view
clip6Display = Show(clip6, renderView1, 'UnstructuredGridRepresentation')

# show color bar/color legend
clip6Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView2.Update()

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(resampleWithDataset1)

# change representation type
resampleWithDataset1Display.SetRepresentationType('Surface With Edges')

# change representation type
resampleWithDataset1Display.SetRepresentationType('Points')

# change representation type
resampleWithDataset1Display.SetRepresentationType('Surface')

# Hide orientation axes
renderView1.OrientationAxesVisibility = 0

# set active view
SetActiveView(lineChartView1)

# set active source
SetActiveSource(plotData1)

# Properties modified on plotData1Display
plotData1Display.UseIndexForXAxis = 1

# set active view
SetActiveView(renderView2)

# set active source
SetActiveSource(resampleWithDataset1)

# set active view
SetActiveView(renderView1)

# set active source
SetActiveSource(clip6)

# hide data in view
Hide(clip6, renderView1)

# set active source
SetActiveSource(clip6)

# show data in view
clip6Display = Show(clip6, renderView1, 'UnstructuredGridRepresentation')

# show color bar/color legend
clip6Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(clip6)

# create a new 'Plane'
plane1 = Plane(registrationName='Plane1')

# show data in view
plane1Display = Show(plane1, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
plane1Display.Representation = 'Surface'
plane1Display.ColorArrayName = [None, '']
plane1Display.SelectTCoordArray = 'TextureCoordinates'
plane1Display.SelectNormalArray = 'Normals'
plane1Display.SelectTangentArray = 'None'
plane1Display.OSPRayScaleArray = 'Normals'
plane1Display.OSPRayScaleFunction = 'PiecewiseFunction'
plane1Display.SelectOrientationVectors = 'None'
plane1Display.ScaleFactor = 0.1
plane1Display.SelectScaleArray = 'None'
plane1Display.GlyphType = 'Arrow'
plane1Display.GlyphTableIndexArray = 'None'
plane1Display.GaussianRadius = 0.005
plane1Display.SetScaleArray = ['POINTS', 'Normals']
plane1Display.ScaleTransferFunction = 'PiecewiseFunction'
plane1Display.OpacityArray = ['POINTS', 'Normals']
plane1Display.OpacityTransferFunction = 'PiecewiseFunction'
plane1Display.DataAxesGrid = 'GridAxesRepresentation'
plane1Display.PolarAxes = 'PolarAxesRepresentation'
plane1Display.SelectInputVectors = ['POINTS', 'Normals']
plane1Display.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
plane1Display.OSPRayScaleFunction.Points = [-0.0155398497892645, 0.0, 0.5, 0.0, -0.0009743210460688503, 0.5, 0.5, 0.0, 0.0135912076971268, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
plane1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 5.878906683738906e-39, 0.5, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
plane1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 5.878906683738906e-39, 0.5, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on plane1
plane1.Point2 = [-0.5, 0.5, 1.0]

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on plane1
plane1.Point2 = [-0.5, 0.0, 1.0]

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on plane1
plane1.Point1 = [1.0, -0.5, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on plane1
plane1.Point2 = [0.0, 0.0, 1.0]

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on plane1
plane1.Point1 = [1.0, 0.0, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on plane1
plane1.Origin = [-0.5, 0.0, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on plane1
plane1.Origin = [0.0, 0.0, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on plane1
plane1.Point1 = [1.0, 0.0, -1.0]

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on plane1
plane1.Origin = [0.0, 0.0, -1.0]

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on plane1
plane1.XResolution = 1000

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on plane1
plane1.YResolution = 300

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on plane1
plane1.XResolution = 300

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Plot On Intersection Curves'
plotOnIntersectionCurves1 = PlotOnIntersectionCurves(registrationName='PlotOnIntersectionCurves1', Input=plane1)
plotOnIntersectionCurves1.SliceType = 'Plane'

# Create a new 'Line Chart View'
lineChartView2 = CreateView('XYChartView')
# uncomment following to set a specific view size
# lineChartView2.ViewSize = [400, 400]

# show data in view
plotOnIntersectionCurves1Display = Show(plotOnIntersectionCurves1, lineChartView2, 'XYChartRepresentation')

# trace defaults for the display properties.
plotOnIntersectionCurves1Display.CompositeDataSetIndex = [0]
plotOnIntersectionCurves1Display.SeriesLabelPrefix = ''

# add view to a layout so it's visible in UI
AssignViewToLayout(view=lineChartView2, layout=layout1, hint=2)

# toggle 3D widget visibility (only when running from the GUI)
Show3DWidgets(proxy=plotOnIntersectionCurves1.SliceType)

# Properties modified on plotOnIntersectionCurves1
plotOnIntersectionCurves1.SliceType = 'Sphere'

# update the view to ensure updated data information
lineChartView2.Update()

# Properties modified on plotOnIntersectionCurves1.SliceType
plotOnIntersectionCurves1.SliceType.Radius = 1.0

# update the view to ensure updated data information
lineChartView2.Update()

# set active source
SetActiveSource(plane1)

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=plotOnIntersectionCurves1.SliceType)

# hide data in view
Hide(plotOnIntersectionCurves1, lineChartView2)

# destroy plotOnIntersectionCurves1
Delete(plotOnIntersectionCurves1)
del plotOnIntersectionCurves1

# update the view to ensure updated data information
lineChartView2.Update()

# destroy lineChartView2
Delete(lineChartView2)
del lineChartView2

# close an empty frame
layout1.Collapse(6)

# set active view
SetActiveView(renderView1)

# set active source
SetActiveSource(resampleWithDataset1)

# set active source
SetActiveSource(plane1)

# Properties modified on plane1
plane1.Point1 = [1.1, 0.0, -1.0]

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on plane1
plane1.Point1 = [1.1, 0.0, -1.1]

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on plane1
plane1.Origin = [0.0, 0.0, -1.01]

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on plane1
plane1.Point1 = [1.1, 0.0, -1.01]

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on plane1
plane1.Point1 = [1.01, 0.0, -1.01]

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on plane1
plane1.Point2 = [0.0, 0.0, 1.01]

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Resample With Dataset'
resampleWithDataset2 = ResampleWithDataset(registrationName='ResampleWithDataset2', SourceDataArrays=clip6,
    DestinationMesh=plane1)
resampleWithDataset2.CellLocator = 'Static Cell Locator'

# show data in view
resampleWithDataset2Display = Show(resampleWithDataset2, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
resampleWithDataset2Display.Representation = 'Surface'
resampleWithDataset2Display.ColorArrayName = [None, '']
resampleWithDataset2Display.SelectTCoordArray = 'None'
resampleWithDataset2Display.SelectNormalArray = 'None'
resampleWithDataset2Display.SelectTangentArray = 'None'
resampleWithDataset2Display.OSPRayScaleArray = 'Analytic'
resampleWithDataset2Display.OSPRayScaleFunction = 'PiecewiseFunction'
resampleWithDataset2Display.SelectOrientationVectors = 'Analytic'
resampleWithDataset2Display.ScaleFactor = 0.20199999809265137
resampleWithDataset2Display.SelectScaleArray = 'Analytic'
resampleWithDataset2Display.GlyphType = 'Arrow'
resampleWithDataset2Display.GlyphTableIndexArray = 'Analytic'
resampleWithDataset2Display.GaussianRadius = 0.010099999904632569
resampleWithDataset2Display.SetScaleArray = ['POINTS', 'Analytic']
resampleWithDataset2Display.ScaleTransferFunction = 'PiecewiseFunction'
resampleWithDataset2Display.OpacityArray = ['POINTS', 'Analytic']
resampleWithDataset2Display.OpacityTransferFunction = 'PiecewiseFunction'
resampleWithDataset2Display.DataAxesGrid = 'GridAxesRepresentation'
resampleWithDataset2Display.PolarAxes = 'PolarAxesRepresentation'
resampleWithDataset2Display.SelectInputVectors = ['POINTS', 'Displacement']
resampleWithDataset2Display.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
resampleWithDataset2Display.OSPRayScaleFunction.Points = [-0.0155398497892645, 0.0, 0.5, 0.0, -0.0009743210460688503, 0.5, 0.5, 0.0, 0.0135912076971268, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
resampleWithDataset2Display.ScaleTransferFunction.Points = [-0.9493356543042809, 0.0, 0.5, 0.0, 0.254239676915044, 0.5, 0.5, 0.0, 1.4578150081343688, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
resampleWithDataset2Display.OpacityTransferFunction.Points = [-0.9493356543042809, 0.0, 0.5, 0.0, 0.254239676915044, 0.5, 0.5, 0.0, 1.4578150081343688, 1.0, 0.5, 0.0]

# hide data in view
Hide(plane1, renderView1)

# hide data in view
Hide(clip6, renderView1)

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(resampleWithDataset2Display, ('POINTS', 'Scalar field'))

# rescale color and/or opacity maps used to include current data range
resampleWithDataset2Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
resampleWithDataset2Display.SetScalarBarVisibility(renderView1, True)

# set active source
SetActiveSource(clip6)

# show data in view
clip6Display = Show(clip6, renderView1, 'UnstructuredGridRepresentation')

# show color bar/color legend
clip6Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# set active view
SetActiveView(lineChartView1)

# destroy lineChartView1
Delete(lineChartView1)
del lineChartView1

RemoveLayout(layout2)

# set active source
SetActiveSource(clip6)

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(resampleWithDataset2)

# set active source
SetActiveSource(plane1)

# set active source
SetActiveSource(plotData1)

# set active source
SetActiveSource(resampleWithDataset1)

# destroy plotData1
Delete(plotData1)
del plotData1

# set active source
SetActiveSource(slice1)

# toggle 3D widget visibility (only when running from the GUI)
Show3DWidgets(proxy=slice1.SliceType)

# hide data in view
Hide(resampleWithDataset1, renderView1)

# show data in view
slice1Display = Show(slice1, renderView1, 'GeometryRepresentation')

# destroy resampleWithDataset1
Delete(resampleWithDataset1)
del resampleWithDataset1

# hide data in view
Hide(resampleWithDataset2, renderView1)

# show data in view
plane1Display = Show(plane1, renderView1, 'GeometryRepresentation')

# destroy resampleWithDataset2
Delete(resampleWithDataset2)
del resampleWithDataset2

# destroy plane1
Delete(plane1)
del plane1

# update the view to ensure updated data information
renderView2.Update()

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(clip9)

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=slice1.SliceType)

# toggle 3D widget visibility (only when running from the GUI)
Show3DWidgets(proxy=clip9.ClipType)

# hide data in view
Hide(slice1, renderView1)

# show data in view
clip9Display = Show(clip9, renderView1, 'UnstructuredGridRepresentation')

# destroy slice1
Delete(slice1)
del slice1

# update the view to ensure updated data information
renderView2.Update()

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(sphere1)

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=clip9.ClipType)

# hide data in view
Hide(clip9, renderView1)

# show data in view
sphere1Display = Show(sphere1, renderView1, 'GeometryRepresentation')

# destroy clip9
Delete(clip9)
del clip9

# update the view to ensure updated data information
renderView2.Update()

# update the view to ensure updated data information
renderView1.Update()

# destroy sphere1
Delete(sphere1)
del sphere1

# update the view to ensure updated data information
renderView1.Update()

# Hide orientation axes
renderView2.OrientationAxesVisibility = 0

# set active view
SetActiveView(renderView1)

# set active source
SetActiveSource(clip8)

# toggle 3D widget visibility (only when running from the GUI)
Show3DWidgets(proxy=clip8.ClipType)

# show data in view
clip8Display = Show(clip8, renderView1, 'UnstructuredGridRepresentation')

# update the view to ensure updated data information
renderView2.Update()

# update the view to ensure updated data information
renderView1.Update()

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=clip8.ClipType)

# reset view to fit data
renderView1.ResetCamera()

#change interaction mode for render view
renderView1.InteractionMode = '2D'

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.75, 4.098051765328873, 0.0]
renderView1.CameraFocalPoint = [0.75, -0.7496136650305125, 0.0]
renderView1.CameraViewUp = [1.0, 0.0, 2.220446049250313e-16]
renderView1.CameraParallelScale = 0.8569552200410602

# current camera placement for renderView2
renderView2.InteractionMode = '2D'
renderView2.CameraPosition = [0.7078228997538776, 9.682673763531604, -0.04963528342434186]
renderView2.CameraFocalPoint = [0.7078228997538776, -0.7498308718204498, -0.04963528342434186]
renderView2.CameraViewUp = [1.0, 0.0, 2.220446049250313e-16]
renderView2.CameraParallelScale = 0.8603449126885949

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).