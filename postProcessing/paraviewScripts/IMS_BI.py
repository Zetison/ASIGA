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
iMS_M4_IE_SSBC_1_meshvtu = XMLUnstructuredGridReader(registrationName='IMS_M4_IE_SSBC_1_mesh.vtu', FileName=['/home/zetison/results/ASIGA/IMS/paraviewResults/IMS_M4_IE_SSBC_1_mesh.vtu'])

# create a new 'XML Unstructured Grid Reader'
iMS_M4_IE_SSBC_2_meshvtu = XMLUnstructuredGridReader(registrationName='IMS_M4_IE_SSBC_2_mesh.vtu', FileName=['/home/zetison/results/ASIGA/IMS/paraviewResults/IMS_M4_IE_SSBC_2_mesh.vtu'])
iMS_M4_IE_SSBC_2_meshvtu.PointArrayStatus = ['Displacement']

# create a new 'XML Unstructured Grid Reader'
iMS_M4_IE_SSBC_1vtu = XMLUnstructuredGridReader(registrationName='IMS_M4_IE_SSBC_1.vtu', FileName=['/home/zetison/results/ASIGA/IMS/paraviewResults/IMS_M4_IE_SSBC_1.vtu'])
iMS_M4_IE_SSBC_1vtu.PointArrayStatus = ['Displacement', 'Scalar field', 'Total scalar field (real)', 'Total scalar field (abs)', 'P_inc']

# create a new 'XML Unstructured Grid Reader'
iMS_M4_IE_SSBC_2vtu = XMLUnstructuredGridReader(registrationName='IMS_M4_IE_SSBC_2.vtu', FileName=['/home/zetison/results/ASIGA/IMS/paraviewResults/IMS_M4_IE_SSBC_2.vtu'])
iMS_M4_IE_SSBC_2vtu.PointArrayStatus = ['Displacement']

# Properties modified on iMS_M4_IE_SSBC_2vtu
iMS_M4_IE_SSBC_2vtu.TimeArray = 'None'

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1920, 1080]

# get layout
layout1 = GetLayout()

# show data in view
iMS_M4_IE_SSBC_2vtuDisplay = Show(iMS_M4_IE_SSBC_2vtu, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
iMS_M4_IE_SSBC_2vtuDisplay.Representation = 'Surface'
iMS_M4_IE_SSBC_2vtuDisplay.ColorArrayName = [None, '']
iMS_M4_IE_SSBC_2vtuDisplay.SelectTCoordArray = 'None'
iMS_M4_IE_SSBC_2vtuDisplay.SelectNormalArray = 'None'
iMS_M4_IE_SSBC_2vtuDisplay.SelectTangentArray = 'None'
iMS_M4_IE_SSBC_2vtuDisplay.OSPRayScaleArray = 'Displacement'
iMS_M4_IE_SSBC_2vtuDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
iMS_M4_IE_SSBC_2vtuDisplay.SelectOrientationVectors = 'Displacement'
iMS_M4_IE_SSBC_2vtuDisplay.ScaleFactor = 7.736840000000003
iMS_M4_IE_SSBC_2vtuDisplay.SelectScaleArray = 'Displacement'
iMS_M4_IE_SSBC_2vtuDisplay.GlyphType = 'Arrow'
iMS_M4_IE_SSBC_2vtuDisplay.GlyphTableIndexArray = 'Displacement'
iMS_M4_IE_SSBC_2vtuDisplay.GaussianRadius = 0.38684200000000013
iMS_M4_IE_SSBC_2vtuDisplay.SetScaleArray = ['POINTS', 'Displacement']
iMS_M4_IE_SSBC_2vtuDisplay.ScaleTransferFunction = 'PiecewiseFunction'
iMS_M4_IE_SSBC_2vtuDisplay.OpacityArray = ['POINTS', 'Displacement']
iMS_M4_IE_SSBC_2vtuDisplay.OpacityTransferFunction = 'PiecewiseFunction'
iMS_M4_IE_SSBC_2vtuDisplay.DataAxesGrid = 'GridAxesRepresentation'
iMS_M4_IE_SSBC_2vtuDisplay.PolarAxes = 'PolarAxesRepresentation'
iMS_M4_IE_SSBC_2vtuDisplay.ScalarOpacityUnitDistance = 2.492080296080013
iMS_M4_IE_SSBC_2vtuDisplay.OpacityArrayName = ['POINTS', 'Displacement']
iMS_M4_IE_SSBC_2vtuDisplay.SelectInputVectors = ['POINTS', 'Displacement']
iMS_M4_IE_SSBC_2vtuDisplay.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
iMS_M4_IE_SSBC_2vtuDisplay.OSPRayScaleFunction.Points = [3.4860574257327247e-09, 0.0, 0.5, 0.0, 2.6920488082154442e-05, 0.0, 0.5, 0.0, 0.00026917350630471283, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
iMS_M4_IE_SSBC_2vtuDisplay.ScaleTransferFunction.Points = [-3.44696467050275e-09, 0.0, 0.5, 0.0, -2.4589021812932597e-09, 0.0, 0.5, 0.0, 6.43366022159215e-09, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
iMS_M4_IE_SSBC_2vtuDisplay.OpacityTransferFunction.Points = [-3.44696467050275e-09, 0.0, 0.5, 0.0, -2.4589021812932597e-09, 0.0, 0.5, 0.0, 6.43366022159215e-09, 1.0, 0.5, 0.0]

# reset view to fit data
renderView1.ResetCamera()

# Properties modified on iMS_M4_IE_SSBC_2_meshvtu
iMS_M4_IE_SSBC_2_meshvtu.TimeArray = 'None'

# show data in view
iMS_M4_IE_SSBC_2_meshvtuDisplay = Show(iMS_M4_IE_SSBC_2_meshvtu, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
iMS_M4_IE_SSBC_2_meshvtuDisplay.Representation = 'Surface'
iMS_M4_IE_SSBC_2_meshvtuDisplay.ColorArrayName = [None, '']
iMS_M4_IE_SSBC_2_meshvtuDisplay.SelectTCoordArray = 'None'
iMS_M4_IE_SSBC_2_meshvtuDisplay.SelectNormalArray = 'None'
iMS_M4_IE_SSBC_2_meshvtuDisplay.SelectTangentArray = 'None'
iMS_M4_IE_SSBC_2_meshvtuDisplay.OSPRayScaleArray = 'Displacement'
iMS_M4_IE_SSBC_2_meshvtuDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
iMS_M4_IE_SSBC_2_meshvtuDisplay.SelectOrientationVectors = 'Displacement'
iMS_M4_IE_SSBC_2_meshvtuDisplay.ScaleFactor = 7.736840000000004
iMS_M4_IE_SSBC_2_meshvtuDisplay.SelectScaleArray = 'Displacement'
iMS_M4_IE_SSBC_2_meshvtuDisplay.GlyphType = 'Arrow'
iMS_M4_IE_SSBC_2_meshvtuDisplay.GlyphTableIndexArray = 'Displacement'
iMS_M4_IE_SSBC_2_meshvtuDisplay.GaussianRadius = 0.3868420000000002
iMS_M4_IE_SSBC_2_meshvtuDisplay.SetScaleArray = ['POINTS', 'Displacement']
iMS_M4_IE_SSBC_2_meshvtuDisplay.ScaleTransferFunction = 'PiecewiseFunction'
iMS_M4_IE_SSBC_2_meshvtuDisplay.OpacityArray = ['POINTS', 'Displacement']
iMS_M4_IE_SSBC_2_meshvtuDisplay.OpacityTransferFunction = 'PiecewiseFunction'
iMS_M4_IE_SSBC_2_meshvtuDisplay.DataAxesGrid = 'GridAxesRepresentation'
iMS_M4_IE_SSBC_2_meshvtuDisplay.PolarAxes = 'PolarAxesRepresentation'
iMS_M4_IE_SSBC_2_meshvtuDisplay.ScalarOpacityUnitDistance = 4.767178108344106
iMS_M4_IE_SSBC_2_meshvtuDisplay.OpacityArrayName = ['POINTS', 'Displacement']
iMS_M4_IE_SSBC_2_meshvtuDisplay.SelectInputVectors = ['POINTS', 'Displacement']
iMS_M4_IE_SSBC_2_meshvtuDisplay.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
iMS_M4_IE_SSBC_2_meshvtuDisplay.OSPRayScaleFunction.Points = [3.4860574257327247e-09, 0.0, 0.5, 0.0, 2.6920488082154442e-05, 0.0, 0.5, 0.0, 0.00026917350630471283, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
iMS_M4_IE_SSBC_2_meshvtuDisplay.ScaleTransferFunction.Points = [-3.4118506259406e-09, 0.0, 0.5, 0.0, -2.427299541187325e-09, 0.0, 0.5, 0.0, 6.4336602215921494e-09, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
iMS_M4_IE_SSBC_2_meshvtuDisplay.OpacityTransferFunction.Points = [-3.4118506259406e-09, 0.0, 0.5, 0.0, -2.427299541187325e-09, 0.0, 0.5, 0.0, 6.4336602215921494e-09, 1.0, 0.5, 0.0]

# Properties modified on iMS_M4_IE_SSBC_1_meshvtu
iMS_M4_IE_SSBC_1_meshvtu.TimeArray = 'None'

# show data in view
iMS_M4_IE_SSBC_1_meshvtuDisplay = Show(iMS_M4_IE_SSBC_1_meshvtu, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
iMS_M4_IE_SSBC_1_meshvtuDisplay.Representation = 'Surface'
iMS_M4_IE_SSBC_1_meshvtuDisplay.ColorArrayName = [None, '']
iMS_M4_IE_SSBC_1_meshvtuDisplay.SelectTCoordArray = 'None'
iMS_M4_IE_SSBC_1_meshvtuDisplay.SelectNormalArray = 'None'
iMS_M4_IE_SSBC_1_meshvtuDisplay.SelectTangentArray = 'None'
iMS_M4_IE_SSBC_1_meshvtuDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
iMS_M4_IE_SSBC_1_meshvtuDisplay.SelectOrientationVectors = 'None'
iMS_M4_IE_SSBC_1_meshvtuDisplay.ScaleFactor = 9.000000000000005
iMS_M4_IE_SSBC_1_meshvtuDisplay.SelectScaleArray = 'None'
iMS_M4_IE_SSBC_1_meshvtuDisplay.GlyphType = 'Arrow'
iMS_M4_IE_SSBC_1_meshvtuDisplay.GlyphTableIndexArray = 'None'
iMS_M4_IE_SSBC_1_meshvtuDisplay.GaussianRadius = 0.45000000000000023
iMS_M4_IE_SSBC_1_meshvtuDisplay.SetScaleArray = [None, '']
iMS_M4_IE_SSBC_1_meshvtuDisplay.ScaleTransferFunction = 'PiecewiseFunction'
iMS_M4_IE_SSBC_1_meshvtuDisplay.OpacityArray = [None, '']
iMS_M4_IE_SSBC_1_meshvtuDisplay.OpacityTransferFunction = 'PiecewiseFunction'
iMS_M4_IE_SSBC_1_meshvtuDisplay.DataAxesGrid = 'GridAxesRepresentation'
iMS_M4_IE_SSBC_1_meshvtuDisplay.PolarAxes = 'PolarAxesRepresentation'
iMS_M4_IE_SSBC_1_meshvtuDisplay.ScalarOpacityUnitDistance = 3.4162504861233285
iMS_M4_IE_SSBC_1_meshvtuDisplay.OpacityArrayName = [None, '']
iMS_M4_IE_SSBC_1_meshvtuDisplay.SelectInputVectors = [None, '']
iMS_M4_IE_SSBC_1_meshvtuDisplay.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
iMS_M4_IE_SSBC_1_meshvtuDisplay.OSPRayScaleFunction.Points = [3.4860574257327247e-09, 0.0, 0.5, 0.0, 2.6920488082154442e-05, 0.0, 0.5, 0.0, 0.00026917350630471283, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
iMS_M4_IE_SSBC_1_meshvtuDisplay.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 0.1, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
iMS_M4_IE_SSBC_1_meshvtuDisplay.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 0.1, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# Properties modified on iMS_M4_IE_SSBC_1vtu
iMS_M4_IE_SSBC_1vtu.TimeArray = 'None'

# show data in view
iMS_M4_IE_SSBC_1vtuDisplay = Show(iMS_M4_IE_SSBC_1vtu, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
iMS_M4_IE_SSBC_1vtuDisplay.Representation = 'Surface'
iMS_M4_IE_SSBC_1vtuDisplay.ColorArrayName = [None, '']
iMS_M4_IE_SSBC_1vtuDisplay.SelectTCoordArray = 'None'
iMS_M4_IE_SSBC_1vtuDisplay.SelectNormalArray = 'None'
iMS_M4_IE_SSBC_1vtuDisplay.SelectTangentArray = 'None'
iMS_M4_IE_SSBC_1vtuDisplay.OSPRayScaleArray = 'Displacement'
iMS_M4_IE_SSBC_1vtuDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
iMS_M4_IE_SSBC_1vtuDisplay.SelectOrientationVectors = 'Displacement'
iMS_M4_IE_SSBC_1vtuDisplay.ScaleFactor = 9.000000000000004
iMS_M4_IE_SSBC_1vtuDisplay.SelectScaleArray = 'Displacement'
iMS_M4_IE_SSBC_1vtuDisplay.GlyphType = 'Arrow'
iMS_M4_IE_SSBC_1vtuDisplay.GlyphTableIndexArray = 'Displacement'
iMS_M4_IE_SSBC_1vtuDisplay.GaussianRadius = 0.4500000000000002
iMS_M4_IE_SSBC_1vtuDisplay.SetScaleArray = ['POINTS', 'Displacement']
iMS_M4_IE_SSBC_1vtuDisplay.ScaleTransferFunction = 'PiecewiseFunction'
iMS_M4_IE_SSBC_1vtuDisplay.OpacityArray = ['POINTS', 'Displacement']
iMS_M4_IE_SSBC_1vtuDisplay.OpacityTransferFunction = 'PiecewiseFunction'
iMS_M4_IE_SSBC_1vtuDisplay.DataAxesGrid = 'GridAxesRepresentation'
iMS_M4_IE_SSBC_1vtuDisplay.PolarAxes = 'PolarAxesRepresentation'
iMS_M4_IE_SSBC_1vtuDisplay.ScalarOpacityUnitDistance = 1.7858721301056923
iMS_M4_IE_SSBC_1vtuDisplay.OpacityArrayName = ['POINTS', 'Displacement']
iMS_M4_IE_SSBC_1vtuDisplay.SelectInputVectors = ['POINTS', 'Displacement']
iMS_M4_IE_SSBC_1vtuDisplay.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
iMS_M4_IE_SSBC_1vtuDisplay.OSPRayScaleFunction.Points = [3.4860574257327247e-09, 0.0, 0.5, 0.0, 2.6920488082154442e-05, 0.0, 0.5, 0.0, 0.00026917350630471283, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
iMS_M4_IE_SSBC_1vtuDisplay.ScaleTransferFunction.Points = [-32808271.4744855, 0.0, 0.5, 0.0, -29527444.327036943, 0.0, 0.5, 0.0, 7.078051567077637e-08, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
iMS_M4_IE_SSBC_1vtuDisplay.OpacityTransferFunction.Points = [-32808271.4744855, 0.0, 0.5, 0.0, -29527444.327036943, 0.0, 0.5, 0.0, 7.078051567077637e-08, 1.0, 0.5, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Clip'
clip1 = Clip(registrationName='Clip1', Input=iMS_M4_IE_SSBC_2vtu)
clip1.ClipType = 'Plane'
clip1.HyperTreeGridClipper = 'Plane'
clip1.Scalars = [None, '']

# init the 'Plane' selected for 'ClipType'
clip1.ClipType.Origin = [-33.99790000000001, 0.0035395015105765815, 0.0]

# init the 'Plane' selected for 'HyperTreeGridClipper'
clip1.HyperTreeGridClipper.Origin = [-33.99790000000001, 0.0035395015105765815, 0.0]

# Properties modified on clip1
clip1.Scalars = ['POINTS', '']

# show data in view
clip1Display = Show(clip1, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
clip1Display.Representation = 'Surface'
clip1Display.ColorArrayName = [None, '']
clip1Display.SelectTCoordArray = 'None'
clip1Display.SelectNormalArray = 'None'
clip1Display.SelectTangentArray = 'None'
clip1Display.OSPRayScaleArray = 'Displacement'
clip1Display.OSPRayScaleFunction = 'PiecewiseFunction'
clip1Display.SelectOrientationVectors = 'Displacement'
clip1Display.ScaleFactor = 3.8684200000000013
clip1Display.SelectScaleArray = 'Displacement'
clip1Display.GlyphType = 'Arrow'
clip1Display.GlyphTableIndexArray = 'Displacement'
clip1Display.GaussianRadius = 0.19342100000000007
clip1Display.SetScaleArray = ['POINTS', 'Displacement']
clip1Display.ScaleTransferFunction = 'PiecewiseFunction'
clip1Display.OpacityArray = ['POINTS', 'Displacement']
clip1Display.OpacityTransferFunction = 'PiecewiseFunction'
clip1Display.DataAxesGrid = 'GridAxesRepresentation'
clip1Display.PolarAxes = 'PolarAxesRepresentation'
clip1Display.ScalarOpacityUnitDistance = 1.6305247581096758
clip1Display.OpacityArrayName = ['POINTS', 'Displacement']
clip1Display.SelectInputVectors = ['POINTS', 'Displacement']
clip1Display.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
clip1Display.OSPRayScaleFunction.Points = [3.4860574257327247e-09, 0.0, 0.5, 0.0, 2.6920488082154442e-05, 0.0, 0.5, 0.0, 0.00026917350630471283, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
clip1Display.ScaleTransferFunction.Points = [-2.62112135616648e-09, 0.0, 0.5, 0.0, -2.1161097797164788e-09, 0.0, 0.5, 0.0, 2.42899440833353e-09, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
clip1Display.OpacityTransferFunction.Points = [-2.62112135616648e-09, 0.0, 0.5, 0.0, -2.1161097797164788e-09, 0.0, 0.5, 0.0, 2.42899440833353e-09, 1.0, 0.5, 0.0]

# hide data in view
Hide(iMS_M4_IE_SSBC_2vtu, renderView1)

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on clip1.ClipType
clip1.ClipType.Normal = [1.0, 1.0, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(iMS_M4_IE_SSBC_2_meshvtu)

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=clip1.ClipType)

# Properties modified on clip1.ClipType
clip1.ClipType.Normal = [0.0, 1.0, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Clip'
clip2 = Clip(registrationName='Clip2', Input=iMS_M4_IE_SSBC_2_meshvtu)
clip2.ClipType = 'Plane'
clip2.HyperTreeGridClipper = 'Plane'
clip2.Scalars = [None, '']

# init the 'Plane' selected for 'ClipType'
clip2.ClipType.Origin = [-33.997900000000016, 0.0035395015105765815, 0.0]

# init the 'Plane' selected for 'HyperTreeGridClipper'
clip2.HyperTreeGridClipper.Origin = [-33.997900000000016, 0.0035395015105765815, 0.0]

# Properties modified on clip2
clip2.Scalars = ['POINTS', '']

# show data in view
clip2Display = Show(clip2, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
clip2Display.Representation = 'Surface'
clip2Display.ColorArrayName = [None, '']
clip2Display.SelectTCoordArray = 'None'
clip2Display.SelectNormalArray = 'None'
clip2Display.SelectTangentArray = 'None'
clip2Display.OSPRayScaleArray = 'Displacement'
clip2Display.OSPRayScaleFunction = 'PiecewiseFunction'
clip2Display.SelectOrientationVectors = 'Displacement'
clip2Display.ScaleFactor = 3.868420000000002
clip2Display.SelectScaleArray = 'Displacement'
clip2Display.GlyphType = 'Arrow'
clip2Display.GlyphTableIndexArray = 'Displacement'
clip2Display.GaussianRadius = 0.1934210000000001
clip2Display.SetScaleArray = ['POINTS', 'Displacement']
clip2Display.ScaleTransferFunction = 'PiecewiseFunction'
clip2Display.OpacityArray = ['POINTS', 'Displacement']
clip2Display.OpacityTransferFunction = 'PiecewiseFunction'
clip2Display.DataAxesGrid = 'GridAxesRepresentation'
clip2Display.PolarAxes = 'PolarAxesRepresentation'
clip2Display.ScalarOpacityUnitDistance = 1.526306133199603
clip2Display.OpacityArrayName = ['POINTS', 'Displacement']
clip2Display.SelectInputVectors = ['POINTS', 'Displacement']
clip2Display.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
clip2Display.OSPRayScaleFunction.Points = [3.4860574257327247e-09, 0.0, 0.5, 0.0, 2.6920488082154442e-05, 0.0, 0.5, 0.0, 0.00026917350630471283, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
clip2Display.ScaleTransferFunction.Points = [-2.62112135616648e-09, 0.0, 0.5, 0.0, -2.1161097797164788e-09, 0.0, 0.5, 0.0, 2.42899440833353e-09, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
clip2Display.OpacityTransferFunction.Points = [-2.62112135616648e-09, 0.0, 0.5, 0.0, -2.1161097797164788e-09, 0.0, 0.5, 0.0, 2.42899440833353e-09, 1.0, 0.5, 0.0]

# hide data in view
Hide(iMS_M4_IE_SSBC_2_meshvtu, renderView1)

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on clip2.ClipType
clip2.ClipType.Normal = [1.0, 1.0, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(iMS_M4_IE_SSBC_1_meshvtu)

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=clip2.ClipType)

# Properties modified on clip2.ClipType
clip2.ClipType.Normal = [0.0, 1.0, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Clip'
clip3 = Clip(registrationName='Clip3', Input=iMS_M4_IE_SSBC_1_meshvtu)
clip3.ClipType = 'Plane'
clip3.HyperTreeGridClipper = 'Plane'
clip3.Scalars = [None, '']

# init the 'Plane' selected for 'ClipType'
clip3.ClipType.Origin = [-33.99790000000001, 0.011329305135951095, 1.7763568394002505e-15]

# init the 'Plane' selected for 'HyperTreeGridClipper'
clip3.HyperTreeGridClipper.Origin = [-33.99790000000001, 0.011329305135951095, 1.7763568394002505e-15]

# Properties modified on clip3
clip3.Scalars = ['POINTS', '']

# show data in view
clip3Display = Show(clip3, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
clip3Display.Representation = 'Surface'
clip3Display.ColorArrayName = [None, '']
clip3Display.SelectTCoordArray = 'None'
clip3Display.SelectNormalArray = 'None'
clip3Display.SelectTangentArray = 'None'
clip3Display.OSPRayScaleFunction = 'PiecewiseFunction'
clip3Display.SelectOrientationVectors = 'None'
clip3Display.ScaleFactor = 4.500000000000003
clip3Display.SelectScaleArray = 'None'
clip3Display.GlyphType = 'Arrow'
clip3Display.GlyphTableIndexArray = 'None'
clip3Display.GaussianRadius = 0.22500000000000012
clip3Display.SetScaleArray = [None, '']
clip3Display.ScaleTransferFunction = 'PiecewiseFunction'
clip3Display.OpacityArray = [None, '']
clip3Display.OpacityTransferFunction = 'PiecewiseFunction'
clip3Display.DataAxesGrid = 'GridAxesRepresentation'
clip3Display.PolarAxes = 'PolarAxesRepresentation'
clip3Display.ScalarOpacityUnitDistance = 1.3029216667591943
clip3Display.OpacityArrayName = [None, '']
clip3Display.SelectInputVectors = [None, '']
clip3Display.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
clip3Display.OSPRayScaleFunction.Points = [3.4860574257327247e-09, 0.0, 0.5, 0.0, 2.6920488082154442e-05, 0.0, 0.5, 0.0, 0.00026917350630471283, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
clip3Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 0.1, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
clip3Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 0.1, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# hide data in view
Hide(iMS_M4_IE_SSBC_1_meshvtu, renderView1)

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on clip3.ClipType
clip3.ClipType.Normal = [1.0, 0.0, 1.0]

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on clip3.ClipType
clip3.ClipType.Normal = [0.0, 0.0, 1.0]

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(iMS_M4_IE_SSBC_1vtu)

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=clip3.ClipType)

# create a new 'Clip'
clip4 = Clip(registrationName='Clip4', Input=iMS_M4_IE_SSBC_1vtu)
clip4.ClipType = 'Plane'
clip4.HyperTreeGridClipper = 'Plane'
clip4.Scalars = ['POINTS', 'P_inc']
clip4.Value = 1.4912156509616636e-08

# init the 'Plane' selected for 'ClipType'
clip4.ClipType.Origin = [-33.997900000000016, 0.011329305135951095, 8.881784197001252e-16]

# init the 'Plane' selected for 'HyperTreeGridClipper'
clip4.HyperTreeGridClipper.Origin = [-33.997900000000016, 0.011329305135951095, 8.881784197001252e-16]

# show data in view
clip4Display = Show(clip4, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
clip4Display.Representation = 'Surface'
clip4Display.ColorArrayName = [None, '']
clip4Display.SelectTCoordArray = 'None'
clip4Display.SelectNormalArray = 'None'
clip4Display.SelectTangentArray = 'None'
clip4Display.OSPRayScaleArray = 'Displacement'
clip4Display.OSPRayScaleFunction = 'PiecewiseFunction'
clip4Display.SelectOrientationVectors = 'Displacement'
clip4Display.ScaleFactor = 4.500000000000002
clip4Display.SelectScaleArray = 'Displacement'
clip4Display.GlyphType = 'Arrow'
clip4Display.GlyphTableIndexArray = 'Displacement'
clip4Display.GaussianRadius = 0.2250000000000001
clip4Display.SetScaleArray = ['POINTS', 'Displacement']
clip4Display.ScaleTransferFunction = 'PiecewiseFunction'
clip4Display.OpacityArray = ['POINTS', 'Displacement']
clip4Display.OpacityTransferFunction = 'PiecewiseFunction'
clip4Display.DataAxesGrid = 'GridAxesRepresentation'
clip4Display.PolarAxes = 'PolarAxesRepresentation'
clip4Display.ScalarOpacityUnitDistance = 1.3904386116265701
clip4Display.OpacityArrayName = ['POINTS', 'Displacement']
clip4Display.SelectInputVectors = ['POINTS', 'Displacement']
clip4Display.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
clip4Display.OSPRayScaleFunction.Points = [3.4860574257327247e-09, 0.0, 0.5, 0.0, 2.6920488082154442e-05, 0.0, 0.5, 0.0, 0.00026917350630471283, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
clip4Display.ScaleTransferFunction.Points = [-32808271.4744855, 0.0, 0.5, 0.0, -29527444.327036943, 0.0, 0.5, 0.0, 7.078051567077637e-08, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
clip4Display.OpacityTransferFunction.Points = [-32808271.4744855, 0.0, 0.5, 0.0, -29527444.327036943, 0.0, 0.5, 0.0, 7.078051567077637e-08, 1.0, 0.5, 0.0]

# hide data in view
Hide(iMS_M4_IE_SSBC_1vtu, renderView1)

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on clip4.ClipType
clip4.ClipType.Normal = [1.0, 0.0, 1.0]

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on clip4.ClipType
clip4.ClipType.Normal = [0.0, 0.0, 1.0]

# update the view to ensure updated data information
renderView1.Update()

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=clip4.ClipType)

# set scalar coloring
ColorBy(clip4Display, ('POINTS', 'Displacement', 'Magnitude'))

# rescale color and/or opacity maps used to include current data range
clip4Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
clip4Display.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'Displacement'
displacementLUT = GetColorTransferFunction('Displacement')

# get opacity transfer function/opacity map for 'Displacement'
displacementPWF = GetOpacityTransferFunction('Displacement')

# set active source
SetActiveSource(clip1)

# toggle 3D widget visibility (only when running from the GUI)
Show3DWidgets(proxy=clip1.ClipType)

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(clip1Display, ('POINTS', 'Displacement', 'Magnitude'))

# rescale color and/or opacity maps used to include current data range
clip1Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
clip1Display.SetScalarBarVisibility(renderView1, True)

# rescale color and/or opacity maps used to exactly fit the current data range
clip1Display.RescaleTransferFunctionToDataRange(False, True)

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=clip1.ClipType)

# hide data in view
Hide(clip3, renderView1)

# hide data in view
Hide(clip2, renderView1)

# set active source
SetActiveSource(clip4)

# Rescale transfer function
displacementLUT.RescaleTransferFunction(1.2993162109553542e-11, 3.69388297967908e+22)

# Rescale transfer function
displacementPWF.RescaleTransferFunction(1.2993162109553542e-11, 3.69388297967908e+22)

# set scalar coloring
ColorBy(clip4Display, ('POINTS', 'Scalar field'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(displacementLUT, renderView1)

# rescale color and/or opacity maps used to include current data range
clip4Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
clip4Display.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'Scalarfield'
scalarfieldLUT = GetColorTransferFunction('Scalarfield')

# get opacity transfer function/opacity map for 'Scalarfield'
scalarfieldPWF = GetOpacityTransferFunction('Scalarfield')

# set active source
SetActiveSource(clip1)

# Rescale transfer function
scalarfieldLUT.RescaleTransferFunction(-3.35900480888283, 1.82124324244669)

# Rescale transfer function
scalarfieldPWF.RescaleTransferFunction(-3.35900480888283, 1.82124324244669)

# rescale color and/or opacity maps used to exactly fit the current data range
clip1Display.RescaleTransferFunctionToDataRange(False, True)

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.CameraPosition = [-29.72536357840762, 37.46930132056051, 96.49801714143355]
renderView1.CameraFocalPoint = [-35.41346161064392, -1.990892948430166, 0.9075438393204439]
renderView1.CameraViewUp = [0.005446533342382111, -0.9244398500158052, 0.3812889966643491]
renderView1.CameraParallelScale = 39.24696798386764

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).