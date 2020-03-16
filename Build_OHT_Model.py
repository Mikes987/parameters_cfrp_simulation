from abaqus import *
from abaqusConstants import *

import numpy as np
import math

import section
import regionToolset
import displayGroupMdbToolset as dgm
import part
import material
import assembly
import step
import interaction
import load
import mesh
import optimization
import job
import sketch
import visualization
import xyPlot
import displayGroupOdbToolset as dgo
import connectorBehavior

def cut(x1,y1,x2,y2,boolean,r, order):

	#--------------------------------------------------------------------------------
	# Function equates linear function made of x1, x2, y1, y2 with circular function
	# made of r.
	# 1. r**2 = x**2 + y**2 <=> y**2 = r**2 - x**2
	# 2. y    = m*x + bl    <=> y**2 = (m*x + bl)**2
	# 3. Equate and use abc equation, then find x, with x find y ==> circle coordinates
	#
	# abc equation contains b**2 - 4ac ==> if its result is >= 0 ==> solution for x
	#--------------------------------------------------------------------------------

	m=(y2-y1)/(x2-x1)
	bl=y1-m*x1

	a=m**2.0+1.0
	b=2.0*m*bl
	c=bl**2.0-r**2.0
	s=b**2.0-4.0*a*c
	if s>=0:
		boolean=False
		lx=(-b+np.sqrt(s))/(2*a)
		ly=m*lx+bl
		return [lx, boolean, ly]
	else:
		boolean=True
		lx=x1
	return [lx, boolean]


#----------------------------------------------------------------------------------
# Functions for partitioning model and look for intersection between partition line
# and circle in the center of the model.
#----------------------------------------------------------------------------------

def upperleft(xh, yh, es, r, a, px1, py1, px2, py2, s):
	px1=px1+es
	dx=xh+px1
	dy=a*dx
	py2=yh-dy
	while py2>-yh:
		s.Line(point1=(px1,py1), point2=(px2,py2))
		results=[px1,py1,px2,py2]
		px1=px1+es
		dx=xh+px1
		dy=a*dx
		py2=yh-dy
	y=yh-0.5*es
	dy=2*yh-0.5*es
	dx=dy/a
	R=2*yh%es
	x=-xh+dx-R
	while y>-yh:
		s.Line(point1=(-xh,y), point2=(x,y))
		y-=0.5*es
		dy-=0.5*es
		dx=dy/a
		x=-xh+dx-R
	return results

def midcheck(xh, yh, es, r, a, boolean, px1, py1, order):
	boolean_check=boolean
	px1=px1+es
	py2=-yh
	dy=2.0*yh
	dx=dy/a
	px2=px1-dx
	while px1<=xh and boolean==boolean_check:
		if order=="circle":
			results=cut(px1,py1,px2,py2,boolean,r, order)
			boolean=results[1]
		else:
			s = order
			if px1 > r + dx / 2 and px1 < dx / 2 + r + 3 * es:
				print 'Partition Line to the right of the hole left out.'
			elif px1 < dx / 2 - r - es and px1 > dx / 2 - r - 3 * es:
				print 'Partition Line to the left of the hole left out.'
			else:
				s.Line(point1=(px1, py1), point2=(px2, py2))
				s.Line(point1=(px1, py1), point2=(px2, py2))
			results=[px1, py1, px2, py2]
		px1+=es
		px2=px1-dx
	return results

def lowerright(xh, yh, es, r, a, px2, s):
	px1=xh
	px2=px2+es
	dx=xh-px2
	dy=a*dx
	py1=-yh+dy
	py2=-yh
	while px2<xh:
		s.Line(point1=(px1,py1), point2=(px2,py2))
		px2+=es
		dx=xh-px2
		dy=a*dx
		py1=-yh+dy
	y=-yh+0.5*es
	dy=2*yh-0.5*es
	dx=dy/a
	x=xh-dx
	while y<yh:
		s.Line(point1=(x,y), point2=(xh,y))
		y+=0.5*es
		dy-=0.5*es
		dx=dy/a
		x=xh-dx

def lowerleft(xh, yh, es, r, a, px1, py1, px2, py2, s):
	px1+=es
	dx=xh+px1
	dy=a*dx
	py2=-yh+dy
	while py2<=yh:
		s.Line(point1=(px1, py1), point2=(-xh, py2))
		results=[px1, py1, px2, py2]
		px1+=es
		dx=xh+px1
		dy=a*dx
		py2=-yh+dy
	y=-yh+0.5*es
	dy=2*yh-0.5*es
	dx=dy/a
	R=2*yh%es
	x=-xh+dx-R
	while y<yh:
		s.Line(point1=(-xh,y), point2=(x, y))
		y+=0.5*es
		dy-=0.5*es
		dx=dy/a
		x=-xh+dx-R
	return results

def midcheck2(xh, yh, es, r, a, boolean, px1, py1, order):
	boolean_check=boolean
	px1=px1+es
	py2=yh
	dy=2*yh
	dx=dy/a
	px2=px1-dx
	while px1<=xh and boolean==boolean_check:
		if order=="circle":
			results=cut(px1,py1,px2,py2,boolean,r, order)
			boolean=results[1]
		else:
			s = order
			if px1 > r + dx / 2 and px1 < dx / 2 + r + 3 * es:
				print 'Partition Line to the right of the hole left out.'
			elif px1 < dx / 2 - r - es and px1 > dx / 2 - r - 3 * es:
				print 'Partition Line to the left of the hole left out'
			else:
				s.Line(point1=(px1, py1), point2=(px2, py2))
			results=[px1, py1, px2, py2]
		px1+=es
		px2=px1-dx
	return results

def upperright(xh, yh, es, r, a, px2, s):
	px1=xh
	px2+=es
	dx=xh-px2
	dy=a*dx
	py1=yh-dy
	py2=yh
	while px2<xh:
		s.Line(point1=(px1, py1), point2=(px2, py2))
		px2+=es
		dx=xh-px2
		dy=a*dx
		py1=yh-dy
	y=yh-0.5*es
	dy=2*yh-0.5*es
	dx=dy/a
	x=xh-dx
	while y>-yh:
		s.Line(point1=(x,y), point2=(xh, y))
		y-=0.5*es
		dy-=0.5*es
		dx=dy/a
		x=xh-dx


#----------------------------------------------------------------------------------
# Main Function "createModel" builds the symmetric model that we need for testing.
#----------------------------------------------------------------------------------

def createModel(angle, numberOfLayer, featureFlags, temperature, es):

	# Begin with an empty sheet, just in Case
	#Mdb()

	# Caution with the upper order. Abaqus crashed twice when it was activated.

	#--------------------------------------------------------------------------------
	# Geometric parameters
	# Only depending, self-typing parameter shall be the orientation angle, given
	# outside of the function.
	#--------------------------------------------------------------------------------

	x_length = 250.0               # entire length of model in x
	y_length =  36.0               # entire length of model in y
	z_length =   0.26              # thickness of each ply

	xh = x_length/2.0              # Needed since we want to set the centerpoint of
	yh = y_length/2.0              # the model at 0,0,0

	hole_diameter = 6.0
	r = hole_diameter/2.0          # needed for modeling and further calculations

	density=1.57e-9
	wkb=0.05

	#--------------------------------------------------------------------------------
	# End of parameters
	# -------------------------------------------------------------------------------

	# Note: element size has to be small for compdam in order to run successfully.

	#--------------------------------------------------------------------------------
	# Creating Part:
	# First locate the model name in abaqus and rename, needed if we use and restart
	# Abaqus multiple times.
	# Plies are created in loops, then assemblied.
	# cf := Carbon Fiber Reinforced Ply
	# -------------------------------------------------------------------------------

	modelname=session.sessionState["Viewport: 1"]["modelName"]
	mdb.models.changeKey(fromName=modelname, toName="OHT")
	cf=mdb.models['OHT']
	#elemenType1 = mesh.ElemType(elemCode=C3D8R, hourglassControl=DEFAULT, distortionControl=DEFAULT,elemLibrary=EXPLICIT)


	#--------------------------------------------------------------------------------
	# Preparing segment for Sections and Materials.
	# 2 Sections:
	# 1. damagable section for the main part of the model.
	# 2. elastic section for the vertical edges of the modes.
	#--------------------------------------------------------------------------------

	vdname="IM7-8552-Damagable"             # VUMAT Damagable name
	vename="IM7-8552-Elastic"               # VUMAT Elastic name
	lname ="Lamina-Elements-"               # Lamina Name

	# Damagable Material
	cf.Material(name=vdname)
	cf.materials[vdname].UserMaterial(mechanicalConstants=( \
			featureFlags, 0, z_length,  0,  0,  0,  0,  0, \
			152689.0, 8703.0,   5164.0,   0.32,     0.45,     62.3,     92.30,    0.277,    \
			0.788,    1.634,    199.8,    0.925,      0,          0,         0,         0,   \
			-5.5e-6,  2.58e-5,  4.06e-9,    5.4,     2326.2,   0.2,      133.3,    0.5, \
			1731.,      0,        0,         0,          0,    wkb,     0,   0.3))
	cf.materials[vdname].Density(table=((density, ), ))

	# Elastic Material
	cf.Material(name=vename)
	cf.materials[vename].UserMaterial(mechanicalConstants=( \
			0, 0, z_length,  0,  0,  0,  0,  0, \
			152689.0, 8703.0,   5164.0,   0.32,     0.45,     62.3,     92.30,    0.277,    \
			0.788,    1.634,    199.8,    0.925,      0,          0,         0,         0,   \
			-5.5e-6,  2.58e-5,  4.06e-9,    5.4,     2326.2,   0.2,      133.3,    0.5, \
			1731.,      0,        0,         0,          0,    wkb,     0,   0.3))
	cf.materials[vename].Density(table=((density, ), ))


	#--------------------------------------------------------------------------------
	# Loop for:
	# - Geometry
	# - Partition
	# - Mesh
	# - Sets
	# - Sections and Properties
	# - Section Assignment
	# - Material Orientation
	#
	#   Included in one Loop, so we can copy a part since it needs to be created more
	#   than ones. Assembly and further creation will appear in the next loop.
	#--------------------------------------------------------------------------------

	i=1
	while i<=numberOfLayer:

		#------------------------------------------------------------------------------
		# First important parameter for all operations: setting a variable for the
		# angle. An angle of 0 means an orientation in x-direction, an angle of 90 an
		# orientation in y-direction.
		#------------------------------------------------------------------------------

		alpha=angle[i-1]
		if alpha!=90.0:
			a=np.tan(np.radians(np.absolute(alpha)))


		#------------------------------------------------------------------------------
		# Do Search here, if a ply with an angle has already been done and copy it. If
		# angles match in numbers, create copy and give it the name "Layer- + (i).
		#------------------------------------------------------------------------------

		if i>1:
			boolean_angle=True
			j=0
			while j<i-1 and boolean_angle==True:
				if alpha==angle[j]:
					name="Layer-"+str(i)
					boolean_angle=False
					p=cf.Part(name=name, objectToCopy=cf.parts["Layer-"+str(j+1)], compressFeatureList=OFF)
					p=cf.parts[name]
					lt="Lamina-Top-"
					lb="Lamina-Bottom-"
					st="Lamina-Elements-"
					print name+" created by copying Layer-"+str(j+1)
					p.surfaces.changeKey(fromName=lt+str(j+1), toName=lt+str(i))
					p.surfaces.changeKey(fromName=lb+str(j+1), toName=lb+str(i))
					p.sets.changeKey(fromName=st+str(j+1), toName=st+str(i))

					#delete old sets
					del p.sets['IM7-8552-Damagable']
					del p.sets['IM7-8552-Elastic']
					#compute new sets
					# Middle ==> damagable
					e_mid = p.elements.getByBoundingBox(-xh + 50, -yh - es, -z_length * 1.1, xh - 50, yh + es,
														z_length * 1.1)
					p.Set(vdname, elements=e_mid)

					# All Elements
					p.Set(elements=p.elements, name=lname + str(i))

					# Edge ==> Make Difference between All Elements and damagable
					p.SetByBoolean(name=vename, operation=DIFFERENCE, sets=(p.sets[lname + str(i)], p.sets[vdname],))

					del p.materialOrientations[0]
					p.MaterialOrientation(region=p.sets[lname + str(i)],
										  orientationType=SYSTEM, axis=AXIS_3, normalAxisDefinition=VECTOR,
										  normalAxisVector=(0.0, 0.0, 1.0), flipNormalDirection=False,
										  normalAxisDirection=AXIS_3, primaryAxisDefinition=VECTOR,
										  primaryAxisVector=(1.0, 0.0, 0.0), primaryAxisDirection=AXIS_1,
										  flipPrimaryDirection=False, additionalRotationType=ROTATION_ANGLE,
										  additionalRotationField='', angle=angle[i - 1], stackDirection=STACK_3)
					session.viewports['Viewport: 1'].setValues(displayedObject=p)
					session.viewports['Viewport: 1'].view.setValues(session.views['Front'])
				else:
					j+=1


		#------------------------------------------------------------------------------
		# 0. Starting to build general part/ply.
		# 1. Create Sketch and draw the rectangle according to x- and y_length
		#------------------------------------------------------------------------------

		if i==1 or boolean_angle==True:
			name="Layer-"+str(i)
			cfp=cf.ConstrainedSketch(name=name ,sheetSize=x_length)
			cfp.rectangle(point1=(-xh,-yh), point2=(xh,yh))


			#----------------------------------------------------------------------------
			# 2. Sketch circle parameter depending on your angle. If the angle is nonzero
			#    and not 90 degrees, find an intersection between a future diagonal
			#    partition line and the circle/hole in order to find proper coordinates
			#    for the hole radius ==> avoiding extra mesh line.
			#    rx, ry := circle radius coordinates
			#----------------------------------------------------------------------------

			if alpha==0.0:
				#es=2.0
				rx=r
				ry=0
			elif alpha==90.0:
				#es=2.0
				rx=0
				ry=r
			elif alpha>0.0:
				#es=2.0
				boolean=True
				py1= yh
				dy=2*yh
				dx=dy/a
				R=dx%es
				px1=-xh+dx+R
				radius=midcheck(xh, yh, es, r, a, boolean, px1, py1, "circle")
				rx=radius[0]
				ry=radius[2]
			else:
				#es=2.0
				boolean=True
				py1=-yh
				dy=2*yh
				dx=dy/a
				R=dx%es
				px1=-xh+dx+R
				radius=midcheck2(xh, yh, es, r, a, boolean, px1, py1, "circle")
				rx=radius[0]
				ry=radius[2]

			cfp.CircleByCenterPerimeter(center=(0.0,0.0), point1=(rx,ry))

			cf.Part(name=name, dimensionality=THREE_D, type=DEFORMABLE_BODY)
			cf.parts[name].BaseSolidExtrude(sketch=cfp, depth=z_length)
			cfp.unsetPrimaryObject()
			p = cf.parts[name]

			session.viewports['Viewport: 1'].setValues(displayedObject=p)
			session.viewports['Viewport: 1'].view.setValues(session.views['Front'])


			#----------------------------------------------------------------------------
			# 3. Create partition sketch for model and do partition according to the
			#    angle. If the angle is 0 or 90 degrees, partitioning can be simply done
			#    by drawing partition lines, for other angles we use specific functions,
			#    which are made for angles <> 0.
			#    - "upperleft":  starting partition at the upper left corner, until the
			#                    lower left corner is reached.
			#    - "midcheck":   continuing with partition until the upper right corner
			#                    is reached.
			#    - "lowerright": continuing and finishing partition until the lower right
			#                    corner is reached.
			#----------------------------------------------------------------------------

			f = p.faces
			t = p.MakeSketchTransform(sketchPlane=f.findAt(coordinates=(xh/2,0.0,z_length)), sketchPlaneSide=SIDE1, origin=(0.0, 0.0, 0.0))
			s = cf.ConstrainedSketch(name='sketch', sheetSize=300.0, gridSpacing=2.5,	transform=t)
			s.setPrimaryObject(option=SUPERIMPOSE)
			p.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)

			if alpha==0.0:
				j=yh
				while j>=-yh:
					s.Line(point1=(-xh, j), point2=(xh, j))
					s.Line(point1=(-xh, j), point2=(xh, j))
					j-=es
				# extra dividing line between circle and the first line
				inter=float(int(r / es)) * es + es / 2.0
				s.Line(point1=(-xh, inter), point2=(xh, inter))
				s.Line(point1 =(-xh, -inter), point2 = (xh, -inter))
			elif alpha==90.0:
				x=-xh
				while x<=xh:
					s.Line(point1=(x,yh), point2=(x,-yh))
					x+=es
				inter=float(int(r / es)) * es + es / 2.0
				# extra dividing line between circle and the first line
				s.Line(point1=(-inter, yh), point2=(-inter, -yh))
				s.Line(point1 =(inter, yh), point2 = (inter, -yh))
			elif alpha>0:
				boolean=True
				px1=-xh
				py1= yh
				px2=-xh
				py2= yh
				coord=upperleft(xh, yh, es, r, a, px1, py1, px2, py2, s)
				px1=coord[0]
				py1=coord[1]
				lx1=coord[0]
				ly1=coord[1]
				lx2=coord[2]
				ly2=coord[3]
				coord=midcheck(xh, yh, es, r, a, boolean, px1, py1, s)
				px2=coord[2]
				rx1=coord[0]
				ry1=coord[1]
				rx2=coord[2]
				ry2=coord[3]
				lowerright(xh, yh, es, r, a, px2, s)
			else:
				boolean=True
				px1=-xh
				py1=-yh
				px2=-xh
				py2=-yh
				coord=lowerleft(xh, yh, es, r, a, px1, py1, px2, py2, s)
				px1=coord[0]
				py1=coord[1]
				lx1=coord[0]
				ly1=coord[1]
				lx2=coord[2]
				ly2=coord[3]
				coord=midcheck2(xh, yh, es, r, a, boolean, px1, py1, s)
				rx1=coord[0]
				ry1=coord[1]
				rx2=coord[2]
				ry2=coord[3]
				px2=coord[2]
				upperright(xh, yh, es, r, a, px2, s)

			pickedFaces = f.findAt(coordinates=(xh/2,0.0,z_length))
			p.PartitionFaceBySketch(faces=pickedFaces, sketch=s)
			s.unsetPrimaryObject()
			del cf.sketches['sketch']


			#----------------------------------------------------------------------------
			# 4. If the angle is not 0 or 90 degrees, we shall split the ply into 3 cells
			#    because we need different mesh-types at the left and right edge and in
			#    the middle of the model in order to not create a distorted/chaotic mesh
			#    at the hole and at the vertical edges.
			#----------------------------------------------------------------------------

			if not (alpha==0.0 or alpha==90.0):
				t = p.MakeSketchTransform(sketchPlane=f.findAt(coordinates=(xh/2,0.0,0.0)), sketchPlaneSide=SIDE2, origin=(0.0, 0.0, 0.0))
				s = cf.ConstrainedSketch(name='sketch', sheetSize=300.0, gridSpacing=2.5,	transform=t)
				s.setPrimaryObject(option=SUPERIMPOSE)
				p.projectReferencesOntoSketch(sketch=s, filter=COPLANAR_EDGES)

				s.Line(point1=(lx1, ly1), point2=(lx2, ly2))
				s.Line(point1=(rx1, ry1), point2=(rx2, ry2))

				pickedFaces = f.findAt(coordinates=(xh/2,0.0,0.0))
				p.PartitionFaceBySketch(faces=pickedFaces, sketch=s)
				s.unsetPrimaryObject()
				del cf.sketches['sketch']
				session.viewports['Viewport: 1'].partDisplay.setValues(mesh=ON)
				session.viewports['Viewport: 1'].partDisplay.meshOptions.setValues(meshTechnique=ON)
				session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(referenceRepresentation=OFF)

			if not (alpha==0.0 or alpha==90.0):
				v=p.vertices
				c=p.cells
				pickedCells=c.findAt(coordinates=(-xh/2, -yh/2, z_length))
				if alpha>0.0:
					p.PartitionCellByPlaneThreePoints(cells=pickedCells, point1=v.findAt(coordinates=(lx1, ly1, 0)), point2=v.findAt(coordinates=(lx1, ly1, z_length)),	point3=v.findAt(coordinates=(lx2, ly2, z_length)))
					pickedCells=c.findAt(coordinates=(xh/2, yh/2, z_length))
					p.PartitionCellByPlaneThreePoints(cells=pickedCells,	point1=v.findAt(coordinates=(rx1, ry1, 0)),	point2=v.findAt(coordinates=(rx1, ry1, z_length)), point3=v.findAt(coordinates=(rx2, ry2, z_length)))
				else:
					p.PartitionCellByPlaneThreePoints(cells=pickedCells, point1=v.findAt(coordinates=(lx1, ly1, 0)), point2=v.findAt(coordinates=(lx1, ly1, z_length)), point3=v.findAt(coordinates=(lx2, ly2, z_length)))
					pickedCells=c.findAt(coordinates=(xh/2, yh/2, z_length))
					p.PartitionCellByPlaneThreePoints(cells=pickedCells, point1=v.findAt(coordinates=(rx1, ry1, 0)), point2=v.findAt(coordinates=(rx1, ry1, z_length)), point3=v.findAt(coordinates=(rx2, ry2, z_length)))


			#----------------------------------------------------------------------------
			# 5. Generate Mesh
			#----------------------------------------------------------------------------

			if alpha==0.0 or alpha==90.0:
				p.seedPart(size=es, minSizeFactor=0.1, deviationFactor=0.6)
				p.setMeshControls(regions=p.cells, elemShape=HEX, algorithm=MEDIAL_AXIS)
				p.generateMesh()
			else:
				c=p.cells
				sin=np.sin(np.arctan(a))
				esmesh=es*sin
				p.seedPart(size=esmesh, minSizeFactor=0.9, deviationFactor=0.8)

				# For angles !=0 or !=90 degrees, further Mesh along the edges is needed to diminish the discrepancy between the
				# element size and mesh size
				Pedges=p.edges.findAt(((0,-y_length/2,0),))
				p.seedEdgeBySize(constraint=FINER, deviationFactor=0.1, edges=Pedges, minSizeFactor=0.9, size=es)
				Pedges = p.edges.findAt(((0, y_length / 2, 0),))
				p.seedEdgeBySize(constraint=FINER, deviationFactor=0.1, edges=Pedges, minSizeFactor=0.9, size=es)


				pickedRegions=c.findAt(coordinates=((0, yh/2, z_length),))
				p.setMeshControls(regions=pickedRegions, elemShape=HEX, algorithm=MEDIAL_AXIS)
				pickedRegions=c.findAt(coordinates=((-xh+es, 0, z_length),))
				p.setMeshControls(regions=pickedRegions, elemShape=HEX_DOMINATED, algorithm=ADVANCING_FRONT)
				pickedRegions=c.findAt(coordinates=((xh-es, 0, z_length),))
				p.setMeshControls(regions=pickedRegions, elemShape=HEX_DOMINATED, algorithm=ADVANCING_FRONT)
				p.generateMesh()


			#----------------------------------------------------------------------------
			# 6. Create Surfaces
			#----------------------------------------------------------------------------

			el=p.elements
			elements=el.getByBoundingBox(-x_length, -y_length, -z_length/2, x_length, y_length, z_length*1.1)
			p.Surface(face2Elements=elements, name="Lamina-Top-"+str(i))
			p.Surface(face1Elements=elements, name="Lamina-Bottom-"+str(i))


			#----------------------------------------------------------------------------
			# 7. Create Sets
			# - At the left and right edge of the model, we have triangles due to the
			#   Partition of the models. We cannot calculate triangles so we have to
			#   "deactivate" that region. Hence, we create a damagable and elastic
			#   material section.
			#----------------------------------------------------------------------------

			# Middle ==> Damagable
			e_mid   = p.elements.getByBoundingBox(-xh + 50, -yh - es, -z_length * 1.1, xh - 50, yh + es, z_length * 1.1)
			p.Set(vdname, elements=e_mid)

			# All Elements
			p.Set(elements=p.elements, name=lname+str(i))

			# Edge ==> Elastic;  Make Difference between All Elements and damagable
			p.SetByBoolean(name=vename, operation=DIFFERENCE, sets=(p.sets[lname+str(i)], p.sets[vdname],))

			# Elementtypes ==>
			#p.setElementType(elemTypes=(mesh.ElemType(elemCode=C3D8R, elemLibrary=EXPLICIT, secondOrderAccuracy=OFF, hourglassControl=ENHANCED), ), regions=(p.sets[lname+str(i)]))
			#p.setElementType(elemTypes=(mesh.ElemType(elemCode=C3D8R, elemLibrary=EXPLICIT, secondOrderAccuracy=OFF, hourglassControl=ENHANCED), ), regions=(p.sets[venmae]))


			#----------------------------------------------------------------------------
			# 8. Creating and Assigning Sections for the Model.
			#----------------------------------------------------------------------------

			cf.HomogeneousSolidSection(name='Section-Damageable', material=vdname, thickness=None)
			p.SectionAssignment(region=p.sets[vdname], sectionName='Section-Damageable', thicknessAssignment=FROM_SECTION)

			cf.HomogeneousSolidSection(name='Section-Elastic', material=vename, thickness=None)
			p.SectionAssignment(region=p.sets[vename], sectionName='Section-Elastic', thicknessAssignment=FROM_SECTION)

			stack_face = p.faces.findAt((hole_diameter, 0.0, z_length), )  # face for stacking ori.
			if alpha == 0.0 or alpha == 90.0:
				p.assignStackDirection(cells=(p.cells[0],), referenceRegion=stack_face)  # assign stacking ori.
			else:
				p.assignStackDirection(cells=(p.cells[0],), referenceRegion=stack_face)  # assign stacking ori.
				p.assignStackDirection(cells=(p.cells[1],), referenceRegion=stack_face)  # assign stacking ori.
				p.assignStackDirection(cells=(p.cells[2],), referenceRegion=stack_face)  # assign stacking ori.


			#----------------------------------------------------------------------------
			# 9. Material Orientation
			#----------------------------------------------------------------------------

			p.MaterialOrientation(region=p.sets[lname+str(i)],
                          orientationType=SYSTEM, axis=AXIS_3, normalAxisDefinition=VECTOR,
                          normalAxisVector=(0.0, 0.0, 1.0), flipNormalDirection=False,
                          normalAxisDirection=AXIS_3, primaryAxisDefinition=VECTOR,
                          primaryAxisVector=(1.0, 0.0, 0.0), primaryAxisDirection=AXIS_1,
                          flipPrimaryDirection=False, additionalRotationType=ROTATION_ANGLE,
                          additionalRotationField='', angle=angle[i - 1], stackDirection=STACK_3)

			print name + " created"

		i += 1


	print "Material Creation done"
	print str(numberOfLayer) + " Parts created"


	#------------------------------------------------------------------------------
	# Create Assembly
	#------------------------------------------------------------------------------

	i=1
	while i<=numberOfLayer:
		name="Layer-"+str(i)
		p=cf.parts[name]
		cfr=cf.rootAssembly
		cfr.DatumCsysByDefault(CARTESIAN)


		#------------------------------------------------------------------------------
		# Setting root Assembly and stacking layers onto each other.
		#------------------------------------------------------------------------------

		cfi=cfr.Instance(name=name, part=p, dependent=ON)

		translate=(0.0, 0.0, z_length*(i-1))
		cfi.translate(translate)

		#------------------------------------------------------------------------------
		# Create Laminate set -hier sollte noch getrennt werden in damagable und non-
		# damagable ansonten ist das gebastel
		#------------------------------------------------------------------------------
		if i == 2:
			cfr.SetByBoolean(name='Laminate', sets=(
				cfr.allInstances['Layer-1'].sets[lname+'1'],
				cfr.allInstances['Layer-2'].sets[lname+'2']))

		if i > 2:
			cfr.SetByBoolean(name='Laminate', sets=(
				cfr.sets['Laminate'],
				cfr.allInstances['Layer-' + str(i)].sets[st + str(i)]))

		if i==2:
			cfr.SetByBoolean(name='Elastic', sets=(
				cfr.allInstances['Layer-1'].sets[vename],
				cfr.allInstances['Layer-2'].sets[vename]))

		if i>2:
			cfr.SetByBoolean(name='Elastic', sets=(
				cfr.sets['Elastic'],
				cfr.allInstances['Layer-'+str(i)].sets[vename]))

		if i==2:
			cfr.SetByBoolean(name='Damagable', sets=(
				cfr.allInstances['Layer-1'].sets[vdname],
				cfr.allInstances['Layer-2'].sets[vdname]))

		if i>2:
			cfr.SetByBoolean(name='Damagable', sets=(
				cfr.sets['Damagable'],
				cfr.allInstances['Layer-'+str(i)].sets[vdname]))

		#-----------------------------------------------------------------------------
		# Nodes on the left edge of the model shall be fixed, nodes on the right edge
		# of the model will be subjected to a load.
		# New Try: Create Fixed-X for x-direction, Fixed-Y for y-direction and Fixed-Z
		# for z-direction
		#-----------------------------------------------------------------------------


		# X-Direction, only Nodes at x = -xh
		nodesfixed_x=cfr.instances[name].nodes.getByBoundingBox(-xh, -yh, z_length*(i-1-0.1), -xh, yh, z_length*i)
		cfr.Set(name="Nodes-Fixed-X-"+str(i), nodes=nodesfixed_x)

		# Y-Direction, only Nodes at y = -yh
		nodesfixed_y=cfr.instances[name].nodes.getByBoundingBox(-xh, -yh, z_length*(i-1-0.1), xh, -yh, z_length*i)
		cfr.Set(name="Nodes-Fixed-Y-"+str(i), nodes=nodesfixed_y)

		# Z-Direction, only Nodes at z = 0
		if i==1:
			nodesfixed_z=cfr.instances[name].nodes.getByBoundingBox(-xh, -yh, z_length*(i-1-0.1), xh, yh, z_length*(i-1+0.1))
			cfr.Set(name="Nodes-Fixed-Z", nodes=nodesfixed_z)

		# Loading Condition
		nodesload=cfr.instances[name].nodes.getByBoundingBox(xh, -yh, z_length*(i-1), xh, yh, z_length*i)
		cfr.Set(name="Nodes-Load-"+str(i), nodes=nodesload)

		#------------------------------------------------------------------------------
		# Summarize all fixed and load nodes
		#------------------------------------------------------------------------------

		if i==2:
			cfr.SetByBoolean(name="Nodes-Fixed-X", sets=(cfr.sets["Nodes-Fixed-X-1"], cfr.sets["Nodes-Fixed-X-2"]))
			cfr.SetByBoolean(name="Nodes-Fixed-Y", sets=(cfr.sets["Nodes-Fixed-Y-1"], cfr.sets["Nodes-Fixed-Y-2"]))
			cfr.SetByBoolean(name="Nodes-Load", sets=(cfr.sets["Nodes-Load-1"], cfr.sets["Nodes-Load-2"]))

		if i>2:
			cfr.SetByBoolean(name="Nodes-Fixed-X", sets=(cfr.sets["Nodes-Fixed-X"], cfr.sets["Nodes-Fixed-X-"+str(i)]))
			cfr.SetByBoolean(name="Nodes-Fixed-Y", sets=(cfr.sets["Nodes-Fixed-Y"], cfr.sets["Nodes-Fixed-Y-"+str(i)]))
			cfr.SetByBoolean(name="Nodes-Load", sets=(cfr.sets["Nodes-Load"], cfr.sets["Nodes-Load-"+str(i)]))
		i+=1


	#--------------------------------------------------------------------------------
	# Reference point is created for easy load application and post processing
	#--------------------------------------------------------------------------------

	cfr = cf.rootAssembly
	rp = cfr.ReferencePoint(point=(x_length*0.8, 0.0, 0.0))
	ref_pt = (cfr.referencePoints.findAt((x_length*0.8, 0.0, 0.0)),)  # find Reference Point in assembly
	cfr.Set(referencePoints=ref_pt, name='RP')
	print "Assembly done"


	#--------------------------------------------------------------------------------
	# Create Equations
	#--------------------------------------------------------------------------------

	cf.Equation(name='RP_Load', terms=((1.0, 'Nodes-Load', 1), (-1.0, 'RP', 1)))

	#--------------------------------------------------------------------------------
	# Create Contact
	#
	# 1. Setting Contact properties
	#--------------------------------------------------------------------------------

	cf.ContactProperty("General")
	cf.interactionProperties["General"].TangentialBehavior(formulation=FRICTIONLESS)
	cf.ContactExp(createStepName="Initial", name="General")
	cf.interactions["General"].includedPairs.setValuesInStep(stepName="Initial", useAllstar=ON)

	cf.ContactProperty("Cohesive-Contact")
	cf.interactionProperties["Cohesive-Contact"].CohesiveBehavior()
	cf.interactionProperties["Cohesive-Contact"].Damage(evolTable=(((0.512, 2.0, 2.0)), ), evolutionType=ENERGY, initTable=((20.0, 80.0, 80.0), ), mixedModeType=BK, useEvolution=ON, useMixedMode=ON)
	cf.interactionProperties["Cohesive-Contact"].damage.setValues(exponent=1.5)

	#--------------------------------------------------------------------------------
	# 2. Connecting the surfaces of the plies.
	#--------------------------------------------------------------------------------

	i=0
	while i<numberOfLayer-1:
		Surface1=cfr.instances["Layer-"+str(i+1)].surfaces["Lamina-Top-"+str(i+1)]
		Surface2=cfr.instances["Layer-"+str(i+2)].surfaces["Lamina-Bottom-"+str(i+2)]
		if i==0:
			cf.interactions["General"].contactPropertyAssignments.appendInStep(assignments=((GLOBAL, SELF, 'General'), (Surface1,Surface2, 'Cohesive-Contact')), stepName='Initial')
		else:
			cf.interactions["General"].contactPropertyAssignments.appendInStep(assignments=((Surface1, Surface2, 'Cohesive-Contact'),), stepName='Initial')
		i+=1

	print "Contacts done"

	#--------------------------------------------------------------------------------
	# Delete single Layer- and nodessets after they have been grouped in the assembly
	#--------------------------------------------------------------------------------

	for j in range(1,numberOfLayer+1):
		del cfr.sets['Nodes-Fixed-X-'+str(j)]
		del cfr.sets['Nodes-Fixed-Y-'+str(j)]
		del cfr.sets['Nodes-Load-'+str(j)]
		print 'Single Set: Nodes-Fixed-X-' +str(j)+ " deleted in Assembly"
		print 'Single Set: Nodes-Fixed-Y-' +str(j)+ " deleted in Assembly"
		print 'Single Set: Nodes-Load-' +str(j)+ " deleted in Assembly"


	#--------------------------------------------------------------------------------
	# Create Boundary Conditions:
	# - Step
	# - Amplitude
	# - Clamp
	# - Load
	#--------------------------------------------------------------------------------

	cf.ExplicitDynamicsStep(name="Step-1", previous="Initial", timePeriod=2.0e-6, massScaling=((SEMI_AUTOMATIC, MODEL, THROUGHOUT_STEP, 0.0, 2.0e-07, BELOW_MIN, 1, 0, 0.0, 0.0, 0, None), ), nlgeom=ON)
	cf.ExplicitDynamicsStep(name="Step-2", previous="Step-1", timePeriod=0.2, massScaling=((SEMI_AUTOMATIC, MODEL, THROUGHOUT_STEP, 0.0, 2.0e-07, BELOW_MIN, 1, 0, 0.0, 0.0, 0, None), ), nlgeom=ON)
	cf.SmoothStepAmplitude(name='Amp-1', timeSpan=STEP, data=((0.0, 0.0), (0.2, 1.0)))
	cf.DisplacementBC(amplitude=UNSET, createStepName='Initial', distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name='Nodes-Fixed-X', region=cfr.sets['Nodes-Fixed-X'], u1=0.0, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)
	cf.DisplacementBC(amplitude=UNSET, createStepName='Initial', distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name='Nodes-Fixed-Y', region=cfr.sets['Nodes-Fixed-Y'], u1=UNSET, u2=0.0, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)
	cf.DisplacementBC(amplitude=UNSET, createStepName='Initial', distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name='Nodes-Fixed-Z', region=cfr.sets['Nodes-Fixed-Z'], u1=UNSET, u2=UNSET, u3=0.0, ur1=UNSET, ur2=UNSET, ur3=UNSET)
	cf.DisplacementBC(amplitude='Amp-1', createStepName='Step-2', distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name='Nodes-Load', region=cfr.sets['RP'], u1=6.0, u2=0.0, u3=0.0, ur1=UNSET, ur2=UNSET, ur3=UNSET)

	print "Boundary Conditions done"


	#--------------------------------------------------------------------------------
	# Set the initital temperature
	# Create a Region. According to the Abaqus macro, it must contain the entire
	# region:
	# All Faces
	# All Cells
	# All edges
	#--------------------------------------------------------------------------------

	edges=[]
	faces=[]
	cells=[]
	vertices=[]

	for i in range(numberOfLayer):
		edges.append(i)
		faces.append(i)
		cells.append(i)
		vertices.append(i)
		edges[i]=cfr.instances["Layer-"+str(i+1)].edges.getByBoundingBox(-xh-1, -yh-1, z_length*i, xh+1, yh+1, z_length*(i+1))
		faces[i]=cfr.instances["Layer-"+str(i+1)].faces.getByBoundingBox(-xh-1, -yh-1, z_length*i, xh+1, yh+1, z_length*(i+1))
		cells[i]=cfr.instances["Layer-"+str(i+1)].cells.getByBoundingBox(-xh-1, -yh-1, z_length*i, xh+1, yh+1, z_length*(i+1))
		vertices[i]=cfr.instances["Layer-"+str(i+1)].vertices.getByBoundingBox(-xh-1, -yh-1, z_length*i, xh+1, yh+1, z_length*(i+1))

	alledges=edges[0]
	allfaces=faces[0]
	allcells=cells[0]
	allvertices=vertices[0]

	if numberOfLayer > 1:
		for i in range (1,numberOfLayer):
			alledges=alledges+edges[i]
			allfaces=allfaces+faces[i]
			allcells=allcells+cells[i]
			allvertices=allvertices+vertices[i]

	region=alledges, allfaces, allcells, allvertices

	#	cfr.Set(alledges, allfaces, allcells, allvertices, name='entire-model')
	# The above statement causes a NameError: name is empty. No Solution found so far.

	cf.Temperature(name='Predefined Field-1', createStepName='Initial', region=region, distributionType=UNIFORM, crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, magnitudes=(temperature,))


	#--------------------------------------------------------------------------------
	# Unsupported Keywords, just taken over by the compdam.py code
	# Doesn't work properly atm hier muessen wir noch ein Element Set hinzufuegen um
	# orphanMeshPart-1.elastic und orphanMeshPart-1.damageable zu ersetzen
    # auf diese Sets koennen wir dann auch die jeweiligen field outputs setzen
	# das heisst sdv usw
	#--------------------------------------------------------------------------------
	cf.keywordBlock.synchVersions(storeNodesAndElements=False)

	position = cf.keywordBlock.sieBlocks.index('*Material, name=IM7-8552-Elastic') + 2
	cf.keywordBlock.insert(position,
			'*Depvar\n' +
			' 26,\n' +
			'  1, CDM_d2\n' +
			'  2, CDM_Fb1\n' +
			'  3, CDM_Fb2\n' +
			'  4, CDM_Fb3\n' +
			'  5, CDM_B\n' +
			'  6, CDM_Lc1\n' +
			'  7, CDM_Lc2\n' +
			'  8, CDM_Lc3\n' +
			'  9, CDM_FIm\n' +
			' 10, CDM_alpha\n' +
			' 11, CDM_STATUS\n' +
			' 12, CDM_Plas12\n' +
			' 13, CDM_Inel12\n' +
			' 14, CDM_FIfT\n' +
			' 15, CDM_slide1\n' +
			' 16, CDM_slide2\n' +
			' 17, CDM_FIfC\n' +
			' 18, CDM_d1T\n' +
			' 19, CDM_d1C\n' +
			' 20, CDM_Plas13\n' +
			' 21, CDM_Inel13\n' +
			' 22, CDM_phi0\n' +
			' 23, CDM_gamma\n' +
			' 24, CDM_Fm1\n' +
			' 25, CDM_Fm2\n' +
			' 26, CDM_Fm3\n' +
			'*Characteristic Length, definition=USER, components=6\n' +
			'*Initial Conditions, Type=Solution\n'+
			' Elastic,  0.d0,  0.d0,  0.d0,  0.d0,  0.d0,  0.d0,  0.d0,\n'+
			' 0.d0,  0.d0,  0,     1,  0.d0,  0.d0,  0.d0,  0.d0,\n'+
			' 0.d0,  0.d0,  0.d0,  0.d0, 0.d0,  0.d0,  0.d0,  0.d0,\n'+
			' 0.d0,  0.d0,  0.d0,')


	position = cf.keywordBlock.sieBlocks.index('*Material, name=IM7-8552-Damagable') + 2
	cf.keywordBlock.insert(position,
			'*Depvar\n' +
			' 26,\n' +
			'  1, CDM_d2\n' +
			'  2, CDM_Fb1\n' +
			'  3, CDM_Fb2\n' +
			'  4, CDM_Fb3\n' +
			'  5, CDM_B\n' +
			'  6, CDM_Lc1\n' +
			'  7, CDM_Lc2\n' +
			'  8, CDM_Lc3\n' +
			'  9, CDM_FIm\n' +
			' 10, CDM_alpha\n' +
			' 11, CDM_STATUS\n' +
			' 12, CDM_Plas12\n' +
			' 13, CDM_Inel12\n' +
			' 14, CDM_FIfT\n' +
			' 15, CDM_slide1\n' +
			' 16, CDM_slide2\n' +
			' 17, CDM_FIfC\n' +
			' 18, CDM_d1T\n' +
			' 19, CDM_d1C\n' +
			' 20, CDM_Plas13\n' +
			' 21, CDM_Inel13\n' +
			' 22, CDM_phi0\n' +
			' 23, CDM_gamma\n' +
			' 24, CDM_Fm1\n' +
			' 25, CDM_Fm2\n' +
			' 26, CDM_Fm3\n' +
			'*Characteristic Length, definition=USER, components=6\n' +
			'*Initial Conditions, Type=Solution\n'+
			' Damagable,  0.d0,  0.d0,  0.d0,  0.d0,  0.d0,  0.d0,  0.d0,\n'+
			' 0.d0,  0.d0,  0,     1,  0.d0,  0.d0,  0.d0,  0.d0,\n'+
			' 0.d0,  0.d0,  0.d0,  0.d0,  1.0d0,  0.d0,  0.d0,  0.d0,\n'+
			' 0.d0,  0.d0,  0.d0,')


def anglesym(angle, featureFlags, temperature, es):

	#--------------------------------------------------------------------------------
	# 1. Flip the given array and save it into a variable
	# 2. Add both arrays to one big symmetric array
	# 3. Save the length of the array in order to know the number of layers
	# 4. Run the symmetric angle and number of layers through the function
	#    "createModel"
	#--------------------------------------------------------------------------------

	rightangle=angle[::-1]
	symangle=angle+rightangle
	l=len(symangle)
	createModel(symangle, l, featureFlags, temperature, es)


#----------------------------------------------------------------------------------
# Starting Point, defining the angles and featureFlags, then the model shall be created by running
# created by running through the functions.
# Parameters:
# - featureFlags: CompDam Parameters
# - angle: first/left set of angles, first function will be "anglesym" where the
#          content of the array will be flipped and saved into a variable. Then,
#          the arrays will be added to create symmetry.
#----------------------------------------------------------------------------------

featureFlags=111001
leftangle = [45.0, -45.0, 0.0, 90.0]
temperature=20
elementsize=1.0

anglesym(leftangle, featureFlags, temperature, elementsize)



#----------------------------------------------------------------------------------
# create field output
#----------------------------------------------------------------------------------
cf = mdb.models['OHT']
cfr = cf.rootAssembly

cf.fieldOutputRequests['F-Output-1'].setValues(
    variables=('S', 'SVAVG', 'PE', 'PEVAVG', 'PEEQ', 'PEEQVAVG', 'LE', 'U', 'V', 'A', 'RF',
               'CSTRESS', 'CSDMG', 'CSMAXSCRT', 'EVF', 'SDV'))

cf.fieldOutputRequests['F-Output-1'].setValues(numIntervals=5)

cf.fieldOutputRequests['F-Output-1'].setValues(
    variables=('S', 'SVAVG', 'PE', 'PEVAVG', 'PEEQ', 'PEEQVAVG', 'LE', 'U', 'V', 'A', 'RF',
               'CSTRESS', 'CSDMG', 'CSMAXSCRT', 'EVF', 'SDV'))

cf.fieldOutputRequests['F-Output-1'].setValues(numIntervals=400)

#----------------------------------------------------------------------------------
# create history output
#----------------------------------------------------------------------------------

regionDef=cfr.sets['RP']
#----------------------------------------------------------------------------------
# safe data
#----------------------------------------------------------------------------------
mdb.Job(model='OHT', name='OHT')
mdb.jobs['OHT'].writeInput()
mdb.saveAs(pathName='./'+'OHT-2019-10-04')
