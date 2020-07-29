# trace generated using paraview version 5.8.0
#
# To ensure correct image size when batch processing, please search 
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'XML Unstructured Grid Reader'
iL_M5_IE_1_meshvtu = XMLUnstructuredGridReader(FileName=['/home/zetison/OneDrive/work/matlab/results/_studies/articleIGA_Ihlenburg3_16/paraviewResults/IL_M5_IE_1_mesh.vtu'])

# create a new 'XML Unstructured Grid Reader'
iL_M5_IE_2_meshvtu = XMLUnstructuredGridReader(FileName=['/home/zetison/OneDrive/work/matlab/results/_studies/articleIGA_Ihlenburg3_16/paraviewResults/IL_M5_IE_2_mesh.vtu'])
iL_M5_IE_2_meshvtu.PointArrayStatus = ['Displacement']

# create a new 'XML Unstructured Grid Reader'
iL_M5_IE_3_meshvtu = XMLUnstructuredGridReader(FileName=['/home/zetison/OneDrive/work/matlab/results/_studies/articleIGA_Ihlenburg3_16/paraviewResults/IL_M5_IE_3_mesh.vtu'])

# create a new 'XML Unstructured Grid Reader'
iL_M5_IE_1vtu = XMLUnstructuredGridReader(FileName=['/home/zetison/OneDrive/work/matlab/results/_studies/articleIGA_Ihlenburg3_16/paraviewResults/IL_M5_IE_1.vtu'])
iL_M5_IE_1vtu.PointArrayStatus = ['Displacement', 'Scalar field', 'Total scalar field (real)', 'Total scalar field (abs)', 'P_inc', 'Analytic', 'Error', 'Error Gradient', 'Error in Energy']

# create a new 'XML Unstructured Grid Reader'
iL_M5_IE_3vtu = XMLUnstructuredGridReader(FileName=['/home/zetison/OneDrive/work/matlab/results/_studies/articleIGA_Ihlenburg3_16/paraviewResults/IL_M5_IE_3.vtu'])
iL_M5_IE_3vtu.PointArrayStatus = ['Displacement', 'Total scalar field (real)', 'Total scalar field (abs)', 'Analytic', 'Error', 'Error Gradient', 'Error in Energy']

# create a new 'XML Unstructured Grid Reader'
iL_M5_IE_2vtu = XMLUnstructuredGridReader(FileName=['/home/zetison/OneDrive/work/matlab/results/_studies/articleIGA_Ihlenburg3_16/paraviewResults/IL_M5_IE_2.vtu'])
iL_M5_IE_2vtu.PointArrayStatus = ['Displacement', 'Analytic', 'Error', 'Error Gradient', 'Error in Energy']

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1216, 751]

# get layout
layout1 = GetLayout()

# show data in view
iL_M5_IE_1vtuDisplay = Show(iL_M5_IE_1vtu, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
iL_M5_IE_1vtuDisplay.Representation = 'Surface'
iL_M5_IE_1vtuDisplay.ColorArrayName = [None, '']
iL_M5_IE_1vtuDisplay.OSPRayScaleArray = 'Analytic'
iL_M5_IE_1vtuDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
iL_M5_IE_1vtuDisplay.SelectOrientationVectors = 'Analytic'
iL_M5_IE_1vtuDisplay.ScaleFactor = 1.2359904729051667
iL_M5_IE_1vtuDisplay.SelectScaleArray = 'Analytic'
iL_M5_IE_1vtuDisplay.GlyphType = 'Arrow'
iL_M5_IE_1vtuDisplay.GlyphTableIndexArray = 'Analytic'
iL_M5_IE_1vtuDisplay.GaussianRadius = 0.061799523645258335
iL_M5_IE_1vtuDisplay.SetScaleArray = ['POINTS', 'Analytic']
iL_M5_IE_1vtuDisplay.ScaleTransferFunction = 'PiecewiseFunction'
iL_M5_IE_1vtuDisplay.OpacityArray = ['POINTS', 'Analytic']
iL_M5_IE_1vtuDisplay.OpacityTransferFunction = 'PiecewiseFunction'
iL_M5_IE_1vtuDisplay.DataAxesGrid = 'GridAxesRepresentation'
iL_M5_IE_1vtuDisplay.PolarAxes = 'PolarAxesRepresentation'
iL_M5_IE_1vtuDisplay.ScalarOpacityUnitDistance = 0.3344997338660675
iL_M5_IE_1vtuDisplay.InputVectors = ['POINTS', 'Displacement']
iL_M5_IE_1vtuDisplay.SelectInputVectors = ['POINTS', 'Displacement']
iL_M5_IE_1vtuDisplay.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
iL_M5_IE_1vtuDisplay.OSPRayScaleFunction.Points = [0.00514701349531455, 0.0, 0.5, 0.0, 0.013980567542770658, 0.0, 0.5, 0.0, 0.0934825539698756, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
iL_M5_IE_1vtuDisplay.ScaleTransferFunction.Points = [-1.20743628534894, 0.0, 0.5, 0.0, -0.869511580758118, 0.0, 0.5, 0.0, 2.17181076055928, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
iL_M5_IE_1vtuDisplay.OpacityTransferFunction.Points = [-1.20743628534894, 0.0, 0.5, 0.0, -0.869511580758118, 0.0, 0.5, 0.0, 2.17181076055928, 1.0, 0.5, 0.0]

# reset view to fit data
renderView1.ResetCamera()

# get the material library
materialLibrary1 = GetMaterialLibrary()

# show data in view
iL_M5_IE_2vtuDisplay = Show(iL_M5_IE_2vtu, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
iL_M5_IE_2vtuDisplay.Representation = 'Surface'
iL_M5_IE_2vtuDisplay.ColorArrayName = [None, '']
iL_M5_IE_2vtuDisplay.OSPRayScaleArray = 'Analytic'
iL_M5_IE_2vtuDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
iL_M5_IE_2vtuDisplay.SelectOrientationVectors = 'Analytic'
iL_M5_IE_2vtuDisplay.ScaleFactor = 1.0149999999999337
iL_M5_IE_2vtuDisplay.SelectScaleArray = 'Analytic'
iL_M5_IE_2vtuDisplay.GlyphType = 'Arrow'
iL_M5_IE_2vtuDisplay.GlyphTableIndexArray = 'Analytic'
iL_M5_IE_2vtuDisplay.GaussianRadius = 0.05074999999999669
iL_M5_IE_2vtuDisplay.SetScaleArray = ['POINTS', 'Analytic']
iL_M5_IE_2vtuDisplay.ScaleTransferFunction = 'PiecewiseFunction'
iL_M5_IE_2vtuDisplay.OpacityArray = ['POINTS', 'Analytic']
iL_M5_IE_2vtuDisplay.OpacityTransferFunction = 'PiecewiseFunction'
iL_M5_IE_2vtuDisplay.DataAxesGrid = 'GridAxesRepresentation'
iL_M5_IE_2vtuDisplay.PolarAxes = 'PolarAxesRepresentation'
iL_M5_IE_2vtuDisplay.ScalarOpacityUnitDistance = 0.3460907782847577
iL_M5_IE_2vtuDisplay.InputVectors = ['POINTS', 'Analytic']
iL_M5_IE_2vtuDisplay.SelectInputVectors = ['POINTS', 'Analytic']
iL_M5_IE_2vtuDisplay.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
iL_M5_IE_2vtuDisplay.OSPRayScaleFunction.Points = [0.00514701349531455, 0.0, 0.5, 0.0, 0.013980567542770658, 0.0, 0.5, 0.0, 0.0934825539698756, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
iL_M5_IE_2vtuDisplay.ScaleTransferFunction.Points = [-8.38674861519538e-11, 0.0, 0.5, 0.0, -6.130067134565662e-11, 0.0, 0.5, 0.0, 1.41800661911018e-10, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
iL_M5_IE_2vtuDisplay.OpacityTransferFunction.Points = [-8.38674861519538e-11, 0.0, 0.5, 0.0, -6.130067134565662e-11, 0.0, 0.5, 0.0, 1.41800661911018e-10, 1.0, 0.5, 0.0]

# show data in view
iL_M5_IE_2_meshvtuDisplay = Show(iL_M5_IE_2_meshvtu, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
iL_M5_IE_2_meshvtuDisplay.Representation = 'Surface'
iL_M5_IE_2_meshvtuDisplay.ColorArrayName = [None, '']
iL_M5_IE_2_meshvtuDisplay.OSPRayScaleArray = 'Displacement'
iL_M5_IE_2_meshvtuDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
iL_M5_IE_2_meshvtuDisplay.SelectOrientationVectors = 'Displacement'
iL_M5_IE_2_meshvtuDisplay.ScaleFactor = 1.0150000000000003
iL_M5_IE_2_meshvtuDisplay.SelectScaleArray = 'Displacement'
iL_M5_IE_2_meshvtuDisplay.GlyphType = 'Arrow'
iL_M5_IE_2_meshvtuDisplay.GlyphTableIndexArray = 'Displacement'
iL_M5_IE_2_meshvtuDisplay.GaussianRadius = 0.05075000000000001
iL_M5_IE_2_meshvtuDisplay.SetScaleArray = ['POINTS', 'Displacement']
iL_M5_IE_2_meshvtuDisplay.ScaleTransferFunction = 'PiecewiseFunction'
iL_M5_IE_2_meshvtuDisplay.OpacityArray = ['POINTS', 'Displacement']
iL_M5_IE_2_meshvtuDisplay.OpacityTransferFunction = 'PiecewiseFunction'
iL_M5_IE_2_meshvtuDisplay.DataAxesGrid = 'GridAxesRepresentation'
iL_M5_IE_2_meshvtuDisplay.PolarAxes = 'PolarAxesRepresentation'
iL_M5_IE_2_meshvtuDisplay.ScalarOpacityUnitDistance = 0.6046762858699317
iL_M5_IE_2_meshvtuDisplay.InputVectors = ['POINTS', 'Displacement']
iL_M5_IE_2_meshvtuDisplay.SelectInputVectors = ['POINTS', 'Displacement']
iL_M5_IE_2_meshvtuDisplay.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
iL_M5_IE_2_meshvtuDisplay.OSPRayScaleFunction.Points = [0.00514701349531455, 0.0, 0.5, 0.0, 0.013980567542770658, 0.0, 0.5, 0.0, 0.0934825539698756, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
iL_M5_IE_2_meshvtuDisplay.ScaleTransferFunction.Points = [-8.54635356272243e-11, 0.0, 0.5, 0.0, -6.262956381829227e-11, 0.0, 0.5, 0.0, 1.4287618246209603e-10, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
iL_M5_IE_2_meshvtuDisplay.OpacityTransferFunction.Points = [-8.54635356272243e-11, 0.0, 0.5, 0.0, -6.262956381829227e-11, 0.0, 0.5, 0.0, 1.4287618246209603e-10, 1.0, 0.5, 0.0]

# show data in view
iL_M5_IE_3vtuDisplay = Show(iL_M5_IE_3vtu, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
iL_M5_IE_3vtuDisplay.Representation = 'Surface'
iL_M5_IE_3vtuDisplay.ColorArrayName = [None, '']
iL_M5_IE_3vtuDisplay.OSPRayScaleArray = 'Analytic'
iL_M5_IE_3vtuDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
iL_M5_IE_3vtuDisplay.SelectOrientationVectors = 'Analytic'
iL_M5_IE_3vtuDisplay.ScaleFactor = 0.9849999999978132
iL_M5_IE_3vtuDisplay.SelectScaleArray = 'Analytic'
iL_M5_IE_3vtuDisplay.GlyphType = 'Arrow'
iL_M5_IE_3vtuDisplay.GlyphTableIndexArray = 'Analytic'
iL_M5_IE_3vtuDisplay.GaussianRadius = 0.04924999999989066
iL_M5_IE_3vtuDisplay.SetScaleArray = ['POINTS', 'Analytic']
iL_M5_IE_3vtuDisplay.ScaleTransferFunction = 'PiecewiseFunction'
iL_M5_IE_3vtuDisplay.OpacityArray = ['POINTS', 'Analytic']
iL_M5_IE_3vtuDisplay.OpacityTransferFunction = 'PiecewiseFunction'
iL_M5_IE_3vtuDisplay.DataAxesGrid = 'GridAxesRepresentation'
iL_M5_IE_3vtuDisplay.PolarAxes = 'PolarAxesRepresentation'
iL_M5_IE_3vtuDisplay.ScalarOpacityUnitDistance = 0.1679307470983999
iL_M5_IE_3vtuDisplay.InputVectors = ['POINTS', 'Displacement']
iL_M5_IE_3vtuDisplay.SelectInputVectors = ['POINTS', 'Displacement']
iL_M5_IE_3vtuDisplay.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
iL_M5_IE_3vtuDisplay.OSPRayScaleFunction.Points = [0.00514701349531455, 0.0, 0.5, 0.0, 0.013980567542770658, 0.0, 0.5, 0.0, 0.0934825539698756, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
iL_M5_IE_3vtuDisplay.ScaleTransferFunction.Points = [-1.91246653161715, 0.0, 0.5, 0.0, -1.401656882009025, 0.0, 0.5, 0.0, 3.1956299644641004, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
iL_M5_IE_3vtuDisplay.OpacityTransferFunction.Points = [-1.91246653161715, 0.0, 0.5, 0.0, -1.401656882009025, 0.0, 0.5, 0.0, 3.1956299644641004, 1.0, 0.5, 0.0]

# show data in view
iL_M5_IE_1_meshvtuDisplay = Show(iL_M5_IE_1_meshvtu, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
iL_M5_IE_1_meshvtuDisplay.Representation = 'Surface'
iL_M5_IE_1_meshvtuDisplay.ColorArrayName = [None, '']
iL_M5_IE_1_meshvtuDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
iL_M5_IE_1_meshvtuDisplay.SelectOrientationVectors = 'None'
iL_M5_IE_1_meshvtuDisplay.ScaleFactor = 1.2359904729056574
iL_M5_IE_1_meshvtuDisplay.SelectScaleArray = 'None'
iL_M5_IE_1_meshvtuDisplay.GlyphType = 'Arrow'
iL_M5_IE_1_meshvtuDisplay.GlyphTableIndexArray = 'None'
iL_M5_IE_1_meshvtuDisplay.GaussianRadius = 0.06179952364528287
iL_M5_IE_1_meshvtuDisplay.SetScaleArray = [None, '']
iL_M5_IE_1_meshvtuDisplay.ScaleTransferFunction = 'PiecewiseFunction'
iL_M5_IE_1_meshvtuDisplay.OpacityArray = [None, '']
iL_M5_IE_1_meshvtuDisplay.OpacityTransferFunction = 'PiecewiseFunction'
iL_M5_IE_1_meshvtuDisplay.DataAxesGrid = 'GridAxesRepresentation'
iL_M5_IE_1_meshvtuDisplay.PolarAxes = 'PolarAxesRepresentation'
iL_M5_IE_1_meshvtuDisplay.ScalarOpacityUnitDistance = 0.5844248659300072
iL_M5_IE_1_meshvtuDisplay.InputVectors = [None, '']
iL_M5_IE_1_meshvtuDisplay.SelectInputVectors = [None, '']
iL_M5_IE_1_meshvtuDisplay.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
iL_M5_IE_1_meshvtuDisplay.OSPRayScaleFunction.Points = [0.00514701349531455, 0.0, 0.5, 0.0, 0.013980567542770658, 0.0, 0.5, 0.0, 0.0934825539698756, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
iL_M5_IE_1_meshvtuDisplay.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 0.1, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
iL_M5_IE_1_meshvtuDisplay.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 0.1, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# show data in view
iL_M5_IE_3_meshvtuDisplay = Show(iL_M5_IE_3_meshvtu, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
iL_M5_IE_3_meshvtuDisplay.Representation = 'Surface'
iL_M5_IE_3_meshvtuDisplay.ColorArrayName = [None, '']
iL_M5_IE_3_meshvtuDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
iL_M5_IE_3_meshvtuDisplay.SelectOrientationVectors = 'None'
iL_M5_IE_3_meshvtuDisplay.ScaleFactor = 0.985
iL_M5_IE_3_meshvtuDisplay.SelectScaleArray = 'None'
iL_M5_IE_3_meshvtuDisplay.GlyphType = 'Arrow'
iL_M5_IE_3_meshvtuDisplay.GlyphTableIndexArray = 'None'
iL_M5_IE_3_meshvtuDisplay.GaussianRadius = 0.04925
iL_M5_IE_3_meshvtuDisplay.SetScaleArray = [None, '']
iL_M5_IE_3_meshvtuDisplay.ScaleTransferFunction = 'PiecewiseFunction'
iL_M5_IE_3_meshvtuDisplay.OpacityArray = [None, '']
iL_M5_IE_3_meshvtuDisplay.OpacityTransferFunction = 'PiecewiseFunction'
iL_M5_IE_3_meshvtuDisplay.DataAxesGrid = 'GridAxesRepresentation'
iL_M5_IE_3_meshvtuDisplay.PolarAxes = 'PolarAxesRepresentation'
iL_M5_IE_3_meshvtuDisplay.ScalarOpacityUnitDistance = 0.2934020401881195
iL_M5_IE_3_meshvtuDisplay.InputVectors = [None, '']
iL_M5_IE_3_meshvtuDisplay.SelectInputVectors = [None, '']
iL_M5_IE_3_meshvtuDisplay.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
iL_M5_IE_3_meshvtuDisplay.OSPRayScaleFunction.Points = [0.00514701349531455, 0.0, 0.5, 0.0, 0.013980567542770658, 0.0, 0.5, 0.0, 0.0934825539698756, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
iL_M5_IE_3_meshvtuDisplay.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 0.1, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
iL_M5_IE_3_meshvtuDisplay.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 0.1, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(iL_M5_IE_1_meshvtu)

# create a new 'Clip'
clip1 = Clip(Input=iL_M5_IE_1_meshvtu)
clip1.ClipType = 'Plane'
clip1.HyperTreeGridClipper = 'Plane'
clip1.Scalars = [None, '']

# init the 'Plane' selected for 'ClipType'
clip1.ClipType.Origin = [0.0, 0.0, 4.440892098500626e-16]

# init the 'Plane' selected for 'HyperTreeGridClipper'
clip1.HyperTreeGridClipper.Origin = [0.0, 0.0, 4.440892098500626e-16]

# Properties modified on clip1
clip1.Scalars = ['POINTS', '']

# show data in view
clip1Display = Show(clip1, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
clip1Display.Representation = 'Surface'
clip1Display.ColorArrayName = [None, '']
clip1Display.OSPRayScaleFunction = 'PiecewiseFunction'
clip1Display.SelectOrientationVectors = 'None'
clip1Display.ScaleFactor = 1.2359904729056574
clip1Display.SelectScaleArray = 'None'
clip1Display.GlyphType = 'Arrow'
clip1Display.GlyphTableIndexArray = 'None'
clip1Display.GaussianRadius = 0.06179952364528287
clip1Display.SetScaleArray = [None, '']
clip1Display.ScaleTransferFunction = 'PiecewiseFunction'
clip1Display.OpacityArray = [None, '']
clip1Display.OpacityTransferFunction = 'PiecewiseFunction'
clip1Display.DataAxesGrid = 'GridAxesRepresentation'
clip1Display.PolarAxes = 'PolarAxesRepresentation'
clip1Display.ScalarOpacityUnitDistance = 0.39988684039562267
clip1Display.InputVectors = [None, '']
clip1Display.SelectInputVectors = [None, '']
clip1Display.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
clip1Display.OSPRayScaleFunction.Points = [0.00514701349531455, 0.0, 0.5, 0.0, 0.013980567542770658, 0.0, 0.5, 0.0, 0.0934825539698756, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
clip1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 0.1, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
clip1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 0.1, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# hide data in view
Hide(iL_M5_IE_1_meshvtu, renderView1)

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on clip1.ClipType
clip1.ClipType.Normal = [1.0, 0.0, 1.0]

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(iL_M5_IE_2_meshvtu)

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=clip1.ClipType)

# Properties modified on clip1.ClipType
clip1.ClipType.Normal = [0.0, 0.0, 1.0]

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Clip'
clip2 = Clip(Input=iL_M5_IE_2_meshvtu)
clip2.ClipType = 'Plane'
clip2.HyperTreeGridClipper = 'Plane'
clip2.Scalars = [None, '']

# init the 'Plane' selected for 'ClipType'
clip2.ClipType.Origin = [0.0, 0.0, 4.440892098500626e-16]

# init the 'Plane' selected for 'HyperTreeGridClipper'
clip2.HyperTreeGridClipper.Origin = [0.0, 0.0, 4.440892098500626e-16]

# Properties modified on clip2
clip2.Scalars = ['POINTS', '']

# show data in view
clip2Display = Show(clip2, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
clip2Display.Representation = 'Surface'
clip2Display.ColorArrayName = [None, '']
clip2Display.OSPRayScaleArray = 'Displacement'
clip2Display.OSPRayScaleFunction = 'PiecewiseFunction'
clip2Display.SelectOrientationVectors = 'Displacement'
clip2Display.ScaleFactor = 1.0150000000000003
clip2Display.SelectScaleArray = 'Displacement'
clip2Display.GlyphType = 'Arrow'
clip2Display.GlyphTableIndexArray = 'Displacement'
clip2Display.GaussianRadius = 0.05075000000000001
clip2Display.SetScaleArray = ['POINTS', 'Displacement']
clip2Display.ScaleTransferFunction = 'PiecewiseFunction'
clip2Display.OpacityArray = ['POINTS', 'Displacement']
clip2Display.OpacityTransferFunction = 'PiecewiseFunction'
clip2Display.DataAxesGrid = 'GridAxesRepresentation'
clip2Display.PolarAxes = 'PolarAxesRepresentation'
clip2Display.ScalarOpacityUnitDistance = 0.4136399310604041
clip2Display.InputVectors = ['POINTS', 'Displacement']
clip2Display.SelectInputVectors = ['POINTS', 'Displacement']
clip2Display.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
clip2Display.OSPRayScaleFunction.Points = [0.00514701349531455, 0.0, 0.5, 0.0, 0.013980567542770658, 0.0, 0.5, 0.0, 0.0934825539698756, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
clip2Display.ScaleTransferFunction.Points = [-8.54635356272243e-11, 0.0, 0.5, 0.0, -6.262956381829227e-11, 0.0, 0.5, 0.0, 1.4287618246209603e-10, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
clip2Display.OpacityTransferFunction.Points = [-8.54635356272243e-11, 0.0, 0.5, 0.0, -6.262956381829227e-11, 0.0, 0.5, 0.0, 1.4287618246209603e-10, 1.0, 0.5, 0.0]

# hide data in view
Hide(iL_M5_IE_2_meshvtu, renderView1)

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on clip2.ClipType
clip2.ClipType.Normal = [1.0, 0.0, 1.0]

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on clip2.ClipType
clip2.ClipType.Normal = [0.0, 0.0, 1.0]

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on clip2.ClipType
clip2.ClipType.Origin = [0.0, 0.0, 1.0]

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(iL_M5_IE_3_meshvtu)

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=clip2.ClipType)

# create a new 'Clip'
clip3 = Clip(Input=iL_M5_IE_3_meshvtu)
clip3.ClipType = 'Plane'
clip3.HyperTreeGridClipper = 'Plane'
clip3.Scalars = [None, '']

# Properties modified on clip3
clip3.Scalars = ['POINTS', '']

# show data in view
clip3Display = Show(clip3, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
clip3Display.Representation = 'Surface'
clip3Display.ColorArrayName = [None, '']
clip3Display.OSPRayScaleFunction = 'PiecewiseFunction'
clip3Display.SelectOrientationVectors = 'None'
clip3Display.ScaleFactor = 0.985
clip3Display.SelectScaleArray = 'None'
clip3Display.GlyphType = 'Arrow'
clip3Display.GlyphTableIndexArray = 'None'
clip3Display.GaussianRadius = 0.04925
clip3Display.SetScaleArray = [None, '']
clip3Display.ScaleTransferFunction = 'PiecewiseFunction'
clip3Display.OpacityArray = [None, '']
clip3Display.OpacityTransferFunction = 'PiecewiseFunction'
clip3Display.DataAxesGrid = 'GridAxesRepresentation'
clip3Display.PolarAxes = 'PolarAxesRepresentation'
clip3Display.ScalarOpacityUnitDistance = 0.20366725142135997
clip3Display.InputVectors = [None, '']
clip3Display.SelectInputVectors = [None, '']
clip3Display.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
clip3Display.OSPRayScaleFunction.Points = [0.00514701349531455, 0.0, 0.5, 0.0, 0.013980567542770658, 0.0, 0.5, 0.0, 0.0934825539698756, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
clip3Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 0.1, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
clip3Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 0.1, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# hide data in view
Hide(iL_M5_IE_3_meshvtu, renderView1)

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
SetActiveSource(iL_M5_IE_1vtu)

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=clip3.ClipType)

# create a new 'Clip'
clip4 = Clip(Input=iL_M5_IE_1vtu)
clip4.ClipType = 'Plane'
clip4.HyperTreeGridClipper = 'Plane'
clip4.Scalars = ['POINTS', 'Analytic']
clip4.Value = 0.4821872376051699

# show data in view
clip4Display = Show(clip4, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
clip4Display.Representation = 'Surface'
clip4Display.ColorArrayName = [None, '']
clip4Display.OSPRayScaleArray = 'Analytic'
clip4Display.OSPRayScaleFunction = 'PiecewiseFunction'
clip4Display.SelectOrientationVectors = 'Analytic'
clip4Display.ScaleFactor = 1.2359904729051667
clip4Display.SelectScaleArray = 'Analytic'
clip4Display.GlyphType = 'Arrow'
clip4Display.GlyphTableIndexArray = 'Analytic'
clip4Display.GaussianRadius = 0.061799523645258335
clip4Display.SetScaleArray = ['POINTS', 'Analytic']
clip4Display.ScaleTransferFunction = 'PiecewiseFunction'
clip4Display.OpacityArray = ['POINTS', 'Analytic']
clip4Display.OpacityTransferFunction = 'PiecewiseFunction'
clip4Display.DataAxesGrid = 'GridAxesRepresentation'
clip4Display.PolarAxes = 'PolarAxesRepresentation'
clip4Display.ScalarOpacityUnitDistance = 0.3649508684864542
clip4Display.InputVectors = ['POINTS', 'Displacement']
clip4Display.SelectInputVectors = ['POINTS', 'Displacement']
clip4Display.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
clip4Display.OSPRayScaleFunction.Points = [0.00514701349531455, 0.0, 0.5, 0.0, 0.013980567542770658, 0.0, 0.5, 0.0, 0.0934825539698756, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
clip4Display.ScaleTransferFunction.Points = [-1.12147598992671, 0.0, 0.5, 0.0, -0.9608014194062812, 0.0, 0.5, 0.0, 0.485269715277578, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
clip4Display.OpacityTransferFunction.Points = [-1.12147598992671, 0.0, 0.5, 0.0, -0.9608014194062812, 0.0, 0.5, 0.0, 0.485269715277578, 1.0, 0.5, 0.0]

# hide data in view
Hide(iL_M5_IE_1vtu, renderView1)

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

# set active source
SetActiveSource(iL_M5_IE_3vtu)

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=clip4.ClipType)

# create a new 'Clip'
clip5 = Clip(Input=iL_M5_IE_3vtu)
clip5.ClipType = 'Plane'
clip5.HyperTreeGridClipper = 'Plane'
clip5.Scalars = ['POINTS', 'Analytic']
clip5.Value = 0.641581716423475

# init the 'Plane' selected for 'ClipType'
clip5.ClipType.Origin = [0.0, 0.0, -8.881784197001252e-16]

# init the 'Plane' selected for 'HyperTreeGridClipper'
clip5.HyperTreeGridClipper.Origin = [0.0, 0.0, -8.881784197001252e-16]

# show data in view
clip5Display = Show(clip5, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
clip5Display.Representation = 'Surface'
clip5Display.ColorArrayName = [None, '']
clip5Display.OSPRayScaleArray = 'Analytic'
clip5Display.OSPRayScaleFunction = 'PiecewiseFunction'
clip5Display.SelectOrientationVectors = 'Analytic'
clip5Display.ScaleFactor = 0.9849999999978132
clip5Display.SelectScaleArray = 'Analytic'
clip5Display.GlyphType = 'Arrow'
clip5Display.GlyphTableIndexArray = 'Analytic'
clip5Display.GaussianRadius = 0.04924999999989066
clip5Display.SetScaleArray = ['POINTS', 'Analytic']
clip5Display.ScaleTransferFunction = 'PiecewiseFunction'
clip5Display.OpacityArray = ['POINTS', 'Analytic']
clip5Display.OpacityTransferFunction = 'PiecewiseFunction'
clip5Display.DataAxesGrid = 'GridAxesRepresentation'
clip5Display.PolarAxes = 'PolarAxesRepresentation'
clip5Display.ScalarOpacityUnitDistance = 0.18315824497382138
clip5Display.InputVectors = ['POINTS', 'Displacement']
clip5Display.SelectInputVectors = ['POINTS', 'Displacement']
clip5Display.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
clip5Display.OSPRayScaleFunction.Points = [0.00514701349531455, 0.0, 0.5, 0.0, 0.013980567542770658, 0.0, 0.5, 0.0, 0.0934825539698756, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
clip5Display.ScaleTransferFunction.Points = [-1.71574873087049, 0.0, 0.5, 0.0, -1.435900823890064, 0.0, 0.5, 0.0, 1.0827303389337697, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
clip5Display.OpacityTransferFunction.Points = [-1.71574873087049, 0.0, 0.5, 0.0, -1.435900823890064, 0.0, 0.5, 0.0, 1.0827303389337697, 1.0, 0.5, 0.0]

# hide data in view
Hide(iL_M5_IE_3vtu, renderView1)

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on clip5.ClipType
clip5.ClipType.Origin = [0.0, 0.0, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on clip5.ClipType
clip5.ClipType.Normal = [1.0, 0.0, 1.0]

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on clip5.ClipType
clip5.ClipType.Normal = [0.0, 0.0, 1.0]

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(iL_M5_IE_2vtu)

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=clip5.ClipType)

# create a new 'Clip'
clip6 = Clip(Input=iL_M5_IE_2vtu)
clip6.ClipType = 'Plane'
clip6.HyperTreeGridClipper = 'Plane'
clip6.Scalars = ['POINTS', 'Error']
clip6.Value = 0.016412016304921808

# init the 'Plane' selected for 'ClipType'
clip6.ClipType.Origin = [0.0, 0.0, -4.440892098500626e-16]

# init the 'Plane' selected for 'HyperTreeGridClipper'
clip6.HyperTreeGridClipper.Origin = [0.0, 0.0, -4.440892098500626e-16]

# show data in view
clip6Display = Show(clip6, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
clip6Display.Representation = 'Surface'
clip6Display.ColorArrayName = [None, '']
clip6Display.OSPRayScaleArray = 'Analytic'
clip6Display.OSPRayScaleFunction = 'PiecewiseFunction'
clip6Display.SelectOrientationVectors = 'Analytic'
clip6Display.ScaleFactor = 1.0149999999999337
clip6Display.SelectScaleArray = 'Analytic'
clip6Display.GlyphType = 'Arrow'
clip6Display.GlyphTableIndexArray = 'Analytic'
clip6Display.GaussianRadius = 0.05074999999999669
clip6Display.SetScaleArray = ['POINTS', 'Analytic']
clip6Display.ScaleTransferFunction = 'PiecewiseFunction'
clip6Display.OpacityArray = ['POINTS', 'Analytic']
clip6Display.OpacityTransferFunction = 'PiecewiseFunction'
clip6Display.DataAxesGrid = 'GridAxesRepresentation'
clip6Display.PolarAxes = 'PolarAxesRepresentation'
clip6Display.ScalarOpacityUnitDistance = 0.3775971019479127
clip6Display.InputVectors = ['POINTS', 'Analytic']
clip6Display.SelectInputVectors = ['POINTS', 'Analytic']
clip6Display.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
clip6Display.OSPRayScaleFunction.Points = [0.00514701349531455, 0.0, 0.5, 0.0, 0.013980567542770658, 0.0, 0.5, 0.0, 0.0934825539698756, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
clip6Display.ScaleTransferFunction.Points = [-8.38674861519538e-11, 0.0, 0.5, 0.0, -6.130067134565662e-11, 0.0, 0.5, 0.0, 1.41800661911018e-10, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
clip6Display.OpacityTransferFunction.Points = [-8.38674861519538e-11, 0.0, 0.5, 0.0, -6.130067134565662e-11, 0.0, 0.5, 0.0, 1.41800661911018e-10, 1.0, 0.5, 0.0]

# hide data in view
Hide(iL_M5_IE_2vtu, renderView1)

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on clip6.ClipType
clip6.ClipType.Origin = [0.0, 0.0, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on clip6.ClipType
clip6.ClipType.Normal = [1.0, 0.0, 1.0]

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on clip6.ClipType
clip6.ClipType.Normal = [0.0, 0.0, 1.0]

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(clip2)

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=clip6.ClipType)

# toggle 3D widget visibility (only when running from the GUI)
Show3DWidgets(proxy=clip2.ClipType)

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(clip6)

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=clip2.ClipType)

# toggle 3D widget visibility (only when running from the GUI)
Show3DWidgets(proxy=clip6.ClipType)

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on clip6.ClipType
clip6.ClipType.Origin = [0.0, 0.0, 1.0]

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(clip6Display, ('POINTS', 'Error in Energy'))

# rescale color and/or opacity maps used to include current data range
clip6Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
clip6Display.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'ErrorinEnergy'
errorinEnergyLUT = GetColorTransferFunction('ErrorinEnergy')

# get opacity transfer function/opacity map for 'ErrorinEnergy'
errorinEnergyPWF = GetOpacityTransferFunction('ErrorinEnergy')

# set active source
SetActiveSource(clip5)

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=clip6.ClipType)

# toggle 3D widget visibility (only when running from the GUI)
Show3DWidgets(proxy=clip5.ClipType)

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(clip5Display, ('POINTS', 'Error in Energy'))

# rescale color and/or opacity maps used to include current data range
clip5Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
clip5Display.SetScalarBarVisibility(renderView1, True)

# set active source
SetActiveSource(clip4)

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=clip5.ClipType)

# toggle 3D widget visibility (only when running from the GUI)
Show3DWidgets(proxy=clip4.ClipType)

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(clip4Display, ('POINTS', 'Error in Energy'))

# rescale color and/or opacity maps used to include current data range
clip4Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
clip4Display.SetScalarBarVisibility(renderView1, True)

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=clip4.ClipType)

# Properties modified on clip4Display
clip4Display.RenderLinesAsTubes = 1

# Properties modified on clip4Display
clip4Display.RenderLinesAsTubes = 0

# set active source
SetActiveSource(iL_M5_IE_3_meshvtu)

# set active source
SetActiveSource(clip3)

# toggle 3D widget visibility (only when running from the GUI)
Show3DWidgets(proxy=clip3.ClipType)

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on clip3Display
clip3Display.RenderLinesAsTubes = 1

# set active source
SetActiveSource(clip2)

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=clip3.ClipType)

# toggle 3D widget visibility (only when running from the GUI)
Show3DWidgets(proxy=clip2.ClipType)

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on clip2Display
clip2Display.RenderLinesAsTubes = 1

# set active source
SetActiveSource(clip1)

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=clip2.ClipType)

# toggle 3D widget visibility (only when running from the GUI)
Show3DWidgets(proxy=clip1.ClipType)

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on clip1Display
clip1Display.RenderLinesAsTubes = 1

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=clip1.ClipType)

# change representation type
clip1Display.SetRepresentationType('Wireframe')

# Hide orientation axes
renderView1.OrientationAxesVisibility = 0

# Properties modified on clip1Display
clip1Display.LineWidth = 2.0

# set active source
SetActiveSource(clip2)

# toggle 3D widget visibility (only when running from the GUI)
Show3DWidgets(proxy=clip2.ClipType)

# update the view to ensure updated data information
renderView1.Update()

# Properties modified on clip2Display
clip2Display.LineWidth = 2.0

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=clip2.ClipType)

# set active source
SetActiveSource(clip3)

# toggle 3D widget visibility (only when running from the GUI)
Show3DWidgets(proxy=clip3.ClipType)

# update the view to ensure updated data information
renderView1.Update()

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=clip3.ClipType)

# Properties modified on clip3Display
clip3Display.LineWidth = 2.0

# change solid color
clip3Display.AmbientColor = [0.0, 0.0, 0.0]
clip3Display.DiffuseColor = [0.0, 0.0, 0.0]

# set active source
SetActiveSource(clip2)

# change solid color
clip2Display.AmbientColor = [0.0, 0.0, 0.0]
clip2Display.DiffuseColor = [0.0, 0.0, 0.0]

# set active source
SetActiveSource(clip1)

# change solid color
clip1Display.AmbientColor = [0.0, 0.0, 0.0]
clip1Display.DiffuseColor = [0.0, 0.0, 0.0]

# Properties modified on clip1Display
clip1Display.LineWidth = 1.0

# set active source
SetActiveSource(clip2)

# Properties modified on clip2Display
clip2Display.LineWidth = 1.0

# set active source
SetActiveSource(clip3)

# Properties modified on clip3Display
clip3Display.LineWidth = 1.0

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.CameraPosition = [0.8948808306178303, -14.614308040420529, 12.71377179871836]
renderView1.CameraFocalPoint = [1.090798890611906, -1.1623927122032074, -1.115242276580657]
renderView1.CameraViewUp = [-0.1116988642971587, 0.7131172195617217, 0.6920890079168173]
renderView1.CameraParallelScale = 10.70399148371416

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).