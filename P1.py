#------------------------------------------------------------------------------
# -*- coding: utf-8 -*-
#------------------------------------------------------------------------------
# -*- coding: utf-8 -*-
"""
P1.py — Abaqus CAE script for Double Shear Bolted Connections

Author:   Md. Ibrahim Kholil (B.Sc. in Civil Engineering)
Email:    engikholil@gmail.com
Version:  0.1.0
Abaqus:   Abaqus 2022
Python:   Abaqus-Python
License:  MIT

Related Work:
"A novel procedure to determine yield and ultimate load and deformation capacity of Double Shear Bolted connections"
(Note: Manuscript under preparation; not yet published)

# DOI: https://doi.org/10.5281/zenodo.16777657

Run with:
abaqus cae noGUI=P1.py

"""


from abaqus import *
from abaqusConstants import *

import __main__
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
import sys
import math

import numpy as np
import matplotlib.pyplot as plt



#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

# functions
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#Instructions
#Input your data in this portions of the code

myJobmodelname = "P1"

myPlateLength = 340.0
myBoltDia = 15.58  #Blot dia
myEndLength =  79.95
myPlateHalfHeight = 50
myPlateThickness =  2.989
myPlateHalfThickness = myPlateThickness/2
myBoltHalfLength = 10.0


myArcX = 0.0      #arc center
myArcY = 0.0
myClearance = 2.12          #Clearance between bolt and plate arc dia
myArcDia = myBoltDia + myClearance   #arc dia
Sh = myBoltHalfLength-myPlateHalfThickness


myDisplacement = 20.0

myDensity = 7.8e-06   #Density
myE = 200000.0  #young modulas
myPoisonRatio = 0.3     #poison ratio
MyFtu = 431.0   #Ultimate strength, Ftu (Mpa)
MyFty = 343.0    #Yield strength, Fty (Mpa)
MySr = 10.7     #Strain at rupture, Ɛr (%)
MyScu = 100      #Ultimate compressive strain
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

def Cae_File_Save(savefilelocation):
    mdb.saveAs(
        pathName=savefilelocation)

#------------------------------------------------------------------------------

def Create_Inp_File(model):
    mdb.jobs[model].writeInput(consistencyChecking=OFF)

#------------------------------------------------------------------------------

def Create_ARC_Point_Solid_line_Create_Plate(model,part,ax,ay,adia,endL,phalfheight,plength,pthickness):
    s1 = mdb.models[model].ConstrainedSketch(name='__profile__', sheetSize=200.0)
    g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
    s1.setPrimaryObject(option=STANDALONE)
    s1.ArcByCenterEnds(center=(ax, ay), point1=(-adia/2, ay), point2=(adia/2, ay),
        direction=CLOCKWISE)
    s1.Spot(point=(-endL, ay))
    s1.Spot(point=(-endL, phalfheight))
    s1.Spot(point=(plength-endL, phalfheight))
    s1.Spot(point=(plength-endL, ay))
    s1.Line(point1=(adia/2, ay), point2=(plength - endL, ay))
    s1.HorizontalConstraint(entity=g[3], addUndoState=False)
    s1.PerpendicularConstraint(entity1=g[2], entity2=g[3], addUndoState=False)
    s1.Line(point1=(plength - endL, ay), point2=(plength - endL, phalfheight))
    s1.VerticalConstraint(entity=g[4], addUndoState=False)
    s1.PerpendicularConstraint(entity1=g[3], entity2=g[4], addUndoState=False)
    s1.Line(point1=(plength - endL, phalfheight), point2=(-endL, phalfheight))
    s1.HorizontalConstraint(entity=g[5], addUndoState=False)
    s1.PerpendicularConstraint(entity1=g[4], entity2=g[5], addUndoState=False)
    s1.Line(point1=(-endL, phalfheight), point2=(-endL, ay))
    s1.VerticalConstraint(entity=g[6], addUndoState=False)
    s1.PerpendicularConstraint(entity1=g[5], entity2=g[6], addUndoState=False)
    s1.Line(point1=(-endL, ay), point2=(-adia/2, ay))
    s1.HorizontalConstraint(entity=g[7], addUndoState=False)
    s1.PerpendicularConstraint(entity1=g[6], entity2=g[7], addUndoState=False)
    p = mdb.models[model].Part(name=part, dimensionality=THREE_D, type=DEFORMABLE_BODY)
    p = mdb.models[model].parts[part]
    p.BaseSolidExtrude(sketch=s1, depth=pthickness)
    s1.unsetPrimaryObject()
    del mdb.models[model].sketches['__profile__']

#------------------------------------------------------------------------------


def Create_Bolt(model,part,ax, ay,bdia,blength):
    s = mdb.models[model].ConstrainedSketch(name='__profile__', 
        sheetSize=200.0)
    g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.setPrimaryObject(option=STANDALONE)
    s.CircleByCenterPerimeter(center=(ax, ay), point1=(-bdia/2, ay))
    p = mdb.models[model].Part(name=part, dimensionality=THREE_D, 
        type=DEFORMABLE_BODY)
    p = mdb.models[model].parts[part]
    p.BaseSolidExtrude(sketch=s, depth=blength)
    s.unsetPrimaryObject()
    del mdb.models[model].sketches['__profile__']

#------------------------------------------------------------------------------


def Plate_material(model,platematerial,density,emodulas,pratio):
    mdb.models[model].Material(name=platematerial)
    mdb.models[model].materials[platematerial].Density(table=((density, ), 
        ))
    mdb.models[model].materials[platematerial].Elastic(table=((emodulas, 
        pratio), ))
    mdb.models[model].materials[platematerial].Plastic(table=(
        (round(True_Stress_8,6), 0.0), 
        (round(True_Stress_9,6),  round(True_Plastic_Strain_9,9)), 
        (round(True_Stress_10,6), round(True_Plastic_Strain_10,9)), 
        (round(True_Stress_11,6), round(True_Plastic_Strain_11,9)), 
        (round(True_Stress_12,6), round(True_Plastic_Strain_12,9)), 
        (round(True_Stress_13,6), round(True_Plastic_Strain_13,9)), 
        (round(True_Stress_14,6), round(True_Plastic_Strain_14,9)), 
        (round(True_Stress_15,6), round(True_Plastic_Strain_15,9)), 
        (round(True_Stress_16,6), round(True_Plastic_Strain_16,9)), 
        (round(True_Stress_17,6), round(True_Plastic_Strain_17,9)), 
        (round(True_Stress_18,6), round(True_Plastic_Strain_18,9)), 
        (round(True_Stress_19,6), round(True_Plastic_Strain_19,9)),
        (round(True_Stress_20,6), round(True_Plastic_Strain_20,9))))
    mdb.models[model].materials[platematerial].DuctileDamageInitiation(
        table=(
            (4.0, -0.33, 0.0), 
            (2.4, 0.05, 0.0), 
            (1.6, 0.1, 0.0), 
            (0.402277608, 0.15, 0.0)))
    #mdb.models[model].materials[platematerial].ductileDamageInitiation.DamageEvolution(
        #type=DISPLACEMENT, table=((1.0, ), ))
    mdb.models[model].materials[platematerial].ShearDamageInitiation(
        ks=-0.2, table=((3.0, 1.65, 0.0), 
        (0.901, 1.731, 0.0), 
        (0.9, 1.732, 0.0), (
        0.901, 1.733, 0.0), 
        (3.0, 1.8, 0.0)))
    #mdb.models[model].materials[platematerial].shearDamageInitiation.DamageEvolution(
        #type=DISPLACEMENT, table=((0.1, ), ))

#------------------------------------------------------------------------------ 


def Bolt_materials(model,boltmaterial,boltdensity,emodulas,pratio):
    mdb.models[model].Material(name=boltmaterial)
    mdb.models[model].materials[boltmaterial].Density(table=((boltdensity, ), 
        ))
    mdb.models[model].materials[boltmaterial].Elastic(table=((emodulas, 
        pratio), ))

#------------------------------------------------------------------------------


def Create_Section(model,crosssection,material):
    mdb.models[model].HomogeneousSolidSection(name=crosssection, material=material, thickness=None)

#------------------------------------------------------------------------------

# Create_Datum_Plane
def Create_Datum_Plane(type_plane,part,model,offset_plane):
    p = mdb.models[model].parts[part]
    myPlane = p.DatumPlaneByPrincipalPlane(principalPlane=type_plane, offset=offset_plane)
    myID = myPlane.id
    return myID

#------------------------------------------------------------------------------


def Create_Partion(model,part,id_plane):
    p = mdb.models[model].parts[part]
    c = p.cells[:]
    d = p.datums
    p.PartitionCellByDatumPlane(datumPlane=d[id_plane], cells=c)

#------------------------------------------------------------------------------

def Section_Assignment(model,part,set_name,section_name):
    p = mdb.models[model].parts[part]
    c = p.cells[:]
    region = p.Set(cells=c, name=set_name)
    p.SectionAssignment(region=region, sectionName=section_name, offset=0.0, 
        offsetType=MIDDLE_SURFACE, offsetField='', 
        thicknessAssignment=FROM_SECTION)

#------------------------------------------------------------------------------

def Assemply(model,part,instance,x,y,z):
    a = mdb.models[model].rootAssembly
    p = mdb.models[model].parts[part]
    a.Instance(name=instance, part=p, dependent=OFF)
    p =a.instances[instance]
    p.translate(vector=(x,y,z))


#------------------------------------------------------------------------------

def Create_Reference_Point(x,y,z,model,setname):
    a = mdb.models[model].rootAssembly
    myRP = a.ReferencePoint(point=(x, y, z))
    r = a.referencePoints
    myRP_Position = r.findAt((x, y, z),)    
    refPoints1=(myRP_Position, )
    a.Set(referencePoints=refPoints1, name=setname)
    return myRP,myRP_Position
#------------------------------------------------------------------------------

def Create_Interaction_Coupling(model,instance,rp_name,bolt_set,rigidbody_name):
    a = mdb.models[model].rootAssembly
    region1=a.sets[rp_name]
    region2=a.instances[instance].sets[bolt_set]
    mdb.models[model].RigidBody(name=rigidbody_name, refPointRegion=region1, bodyRegion=region2)
#------------------------------------------------------------------------------



def Create_Step(model,pre_step_name,step_name,nlg_on_off,max_num_inc,initial_inc,min_inc,max_inc):
    a = mdb.models[model].ImplicitDynamicsStep(name=step_name, previous=pre_step_name, application=QUASI_STATIC, nohaf=OFF, amplitude=RAMP, alpha=DEFAULT, initialConditions=OFF, nlgeom=nlg_on_off)
    a = mdb.models[model].steps[step_name].setValues(maxNumInc=max_num_inc, initialInc=initial_inc, minInc=min_inc, maxInc=max_inc)

#------------------------------------------------------------------------------


def Fieldoutput_History_Output_Request(model,fieldoutput,historyoutput,timeintervel):
    a = mdb.models[model].fieldOutputRequests[fieldoutput].setValues(variables=('S', 'PE', 'PEEQ', 'PEMAG', 'LE', 'U', 'V', 'A', 'RF', 'CF', 'CSTRESS', 'CDISP', 'EVOL'), timeInterval=timeintervel, timeMarks=OFF)
    a = mdb.models[model].historyOutputRequests[historyoutput].setValues(timeInterval=timeintervel, timeMarks=OFF)

#------------------------------------------------------------------------------


def Contact_Property(model,contactprop,frictionfactor):
    mdb.models[model].ContactProperty(contactprop)
    mdb.models[model].interactionProperties[contactprop].NormalBehavior(pressureOverclosure=HARD, allowSeparation=ON, constraintEnforcementMethod=DEFAULT)
    mdb.models[model].interactionProperties[contactprop].TangentialBehavior(formulation=PENALTY, directionality=ISOTROPIC, slipRateDependency=OFF, pressureDependency=OFF, temperatureDependency=OFF, dependencies=0, table=((frictionfactor, ), ), shearStressLimit=None, maximumElasticSlip=FRACTION, fraction=0.005, elasticSlipStiffness=None)

#------------------------------------------------------------------------------


def Bolt_Plate_Surface(model,instance1,instance2,boltsurface1,platesurface1):
    a = mdb.models[model].rootAssembly
    s1 = a.instances[instance1].faces
    side1Faces1 = s1.getSequenceFromMask(mask=('[#3014 ]', ), )
    a.Surface(side1Faces=side1Faces1, name=boltsurface1)
    s2 = a.instances[instance2].faces
    side2Faces1 = s2.getSequenceFromMask(mask=('[#210 ]', ), )
    a.Surface(side1Faces=side2Faces1, name=platesurface1)

#------------------------------------------------------------------------------



def Create_Contact_BP_Surface(model,boltsurface1,platesurface1,bpcontact,step,contactprop):
    a = mdb.models[model].rootAssembly
    region1=a.surfaces[boltsurface1]
    region2=a.surfaces[platesurface1]
    mdb.models[model].SurfaceToSurfaceContactStd(name=bpcontact,createStepName=step, main=region1, secondary=region2, sliding=FINITE,thickness=ON, interactionProperty=contactprop, adjustMethod=OVERCLOSED,initialClearance=OMIT, datumAxis=None, clearanceRegion=None, tied=OFF)


#------------------------------------------------------------------------------


def Bolt_Displacement_Boundary_Conditions(model,set_name,bc_name,step_name,ux,uy,uz,urx,ury,urz):
    a = mdb.models[model].rootAssembly
    region = a.sets[set_name]
    mdb.models[model].DisplacementBC(name=bc_name, createStepName=step_name, region=region, u1=ux, u2=uy, u3=uz, ur1=urx, ur2=ury, ur3=urz, amplitude=UNSET, fixed=OFF,  distributionType=UNIFORM, fieldName='', localCsys=None)

#------------------------------------------------------------------------------


def Create_Surface(model,part,x,y,z,surface_name):
    a = mdb.models[model].rootAssembly
    s1 = a.instances[part].faces
    side1Faces1 = s1.findAt(((x,y,z), ))
    a.Surface(side1Faces=side1Faces1, name=surface_name)

#------------------------------------------------------------------------------


def Create_BC_Plate_Edge_Fixed(model,instance,plateedge,bc_name,step_name):
    a = mdb.models[model].rootAssembly
    f1 = a.instances[instance].faces
    faces1 = f1.getSequenceFromMask(mask=('[#80 ]', ), )
    region = a.Set(faces=faces1, name=plateedge)
    mdb.models[model].EncastreBC(name=bc_name, 
        createStepName=step_name, region=region, localCsys=None)

#------------------------------------------------------------------------------

def Create_Extra_BC_For_Thin_Plate_Edge_Zaxis_Fixed(model,instance,surface_set_name,bc_name,x,y,z):
    a = mdb.models[model].rootAssembly
    f1 = a.instances[instance].faces
    faces1 = f1.findAt(((x,y,z), ))
    region = a.Set(faces=faces1, name=surface_set_name)
    mdb.models[model].DisplacementBC(name=bc_name,createStepName='Loading', region=region, u1=UNSET, u2=UNSET, u3=0.0,ur1=UNSET, ur2=UNSET, ur3=UNSET, amplitude=UNSET, fixed=OFF,distributionType=UNIFORM, fieldName='', localCsys=None)

#------------------------------------------------------------------------------

def Create_Symmetry(model,instance,platebottom,sym,step_name):
    a = mdb.models[model].rootAssembly
    f1 = a.instances[instance].faces
    faces1 = f1.getSequenceFromMask(mask=('[#500 ]', ), )
    region = a.Set(faces=faces1, name=platebottom)
    mdb.models[model].YsymmBC(name=sym, createStepName=step_name, region=region, localCsys=None)

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
def Create_Z_Symmetry(model,instance_1,instance_2,set_name,sym_name):
    a = mdb.models[model].rootAssembly
    f1 = a.instances[instance_1].faces
    faces1 = f1.getSequenceFromMask(mask=('[#1002 ]', ), )
    f2 = a.instances[instance_2].faces
    faces2 = f2.getSequenceFromMask(mask=('[#8460 ]', ), )
    region = a.Set(faces=faces1+faces2, name=set_name)
    mdb.models[model].ZsymmBC(name=sym_name, 
        createStepName='Loading', region=region, localCsys=None)
    

#------------------------------------------------------------------------------


#Create Mesh
def Create_Mesh(model,instance_1,instance_2,mesh_size,edge_size_1):
    a = mdb.models[model].rootAssembly
    partInstances =(a.instances[instance_1], a.instances[instance_2], )
    e1 = a.instances[instance_1].edges
    pickedEdges = e1.getSequenceFromMask(mask=('[#14050 ]', ), )
    a.seedEdgeBySize(edges=pickedEdges, size=edge_size_1, deviationFactor=0.1, 
        minSizeFactor=0.1, constraint=FINER)
    a.seedPartInstance(regions=partInstances, size=mesh_size, deviationFactor=0.1, minSizeFactor=0.1)
    a.generateMesh(regions=partInstances)


#------------------------------------------------------------------------------


def Select_Element_Type(model,instance_1,instance_2,element_del,maxdeg):
    elemType1 = mesh.ElemType(elemCode=C3D8R, elemLibrary=STANDARD, kinematicSplit=AVERAGE_STRAIN, secondOrderAccuracy=OFF, hourglassControl=DEFAULT, distortionControl=DEFAULT, elemDeletion=element_del, maxDegradation=maxdeg)
    elemType2 = mesh.ElemType(elemCode=C3D6, elemLibrary=STANDARD)
    elemType3 = mesh.ElemType(elemCode=C3D4, elemLibrary=STANDARD)
    a = mdb.models[model].rootAssembly
    c1 = a.instances[instance_1].cells
    c2 = a.instances[instance_2].cells
    pickedRegions =((c1+c2), )
    a.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2, elemType3))


#------------------------------------------------------------------------------

def Create_Node_Set_ByBoundingBox(model,instance,x1,y1,z1,x2,y2,z2,set_name):
    a = mdb.models[model].rootAssembly
    n1 = a.instances[instance].nodes
    nodes1 = n1.getByBoundingBox(x1,y1,z1,x2,y2,z2)
    a.Set(nodes=nodes1, name=set_name)

      
#----------------------------------------------------------------------------

def Create_Element_Set_ByBoundingBox(model,instance,x1,y1,z1,x2,y2,z2,set_name):
    a = mdb.models[model].rootAssembly
    e1 = a.instances[instance].elements
    elements1 = e1.getByBoundingBox(x1,y1,z1,x2,y2,z2)
    a.Set(elements=elements1, name=set_name)

#----------------------------------------------------------------------------

def Create_Job(model,job_name, cpu):
    a = mdb.models[model].rootAssembly
    mdb.Job(name=job_name, model=model, description='', type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=cpu, numGPUs=0)

#------------------------------------------------------------------------------

#Job submit
def Job_Submit(myJobName):
    os.chdir(myDirectory)
    mdb.jobs[myJobName].submit(consistencyChecking=OFF)
    

#------------------------------------------------------------------------------ 



myPart_1 = "Plate"
myPart_2 = "Bolt"
myString = myJobmodelname
mdb.Model(name=myString)
myMaterialName_1 = "Plate"
myMaterialName_2 = "Bolt"
myPlateCs = "Plate CS"
myBoltCs = "Bolt CS"
myJobName= myJobmodelname



#------------------------------------------------------------------------------ 

#myfilesave = 'H:/Quater model By Aziz Sir/PS_4/P1 '
#myDirectory = r"H:\Quater model By Aziz Sir\PS_9"
#------------------------------------------------------------------------------    

MySus = 100.0*((MySr/100.0)-(MyFtu/myE))   #Uniform strain at rupture, Ɛus (%)
Mysus_1 = MySus/0.2
MyMath_1 = (np.log(Mysus_1))
#print(MyMath_1)
Ratioft = float(MyFtu/MyFty)
#print(Ratioft)
MyMath_2 = (np.log(Ratioft))  #Strain hardening exponent, n
#print(MyMath_2)
Myn = (MyMath_1/MyMath_2)
#print(Myn)

i=0

#Engg_Stress
Engg_Stress_1 = i*MyFty/5
Engg_Stress_2 = Engg_Stress_1 + 1*MyFty/5
Engg_Stress_3 = Engg_Stress_1 + 2*MyFty/5
Engg_Stress_4 = Engg_Stress_1 + 3*MyFty/5
Engg_Stress_5 = Engg_Stress_1 + 4*MyFty/5


#Engg_Strain
Engg_Strain_1 = Engg_Stress_1/myE+0.002*(Engg_Stress_1/MyFty)**Myn
Engg_Strain_2 = Engg_Stress_2/myE+0.002*(Engg_Stress_2/MyFty)**Myn
Engg_Strain_3 = Engg_Stress_3/myE+0.002*(Engg_Stress_3/MyFty)**Myn
Engg_Strain_4 = Engg_Stress_4/myE+0.002*(Engg_Stress_4/MyFty)**Myn
Engg_Strain_5 = Engg_Stress_5/myE+0.002*(Engg_Stress_5/MyFty)**Myn

#True_Stress
True_Stress_1 = Engg_Stress_1*(1+Engg_Strain_1)
True_Stress_2 = Engg_Stress_2*(1+Engg_Strain_2)
True_Stress_3 = Engg_Stress_3*(1+Engg_Strain_3)
True_Stress_4 = Engg_Stress_4*(1+Engg_Strain_4)
True_Stress_5 = Engg_Stress_5*(1+Engg_Strain_5)

#True_Plastic_Strain
True_Plastic_Strain_1 = np.log(1+Engg_Strain_1) - True_Stress_1/myE
True_Plastic_Strain_2 = np.log(1+Engg_Strain_2) - True_Stress_2/myE
True_Plastic_Strain_3 = np.log(1+Engg_Strain_3) - True_Stress_3/myE
True_Plastic_Strain_4 = np.log(1+Engg_Strain_4) - True_Stress_4/myE
True_Plastic_Strain_5 = np.log(1+Engg_Strain_5) - True_Stress_5/myE


#print(round(True_Plastic_Strain_1,6),round(True_Plastic_Strain_2,6),round(True_Plastic_Strain_3,6),round(True_Plastic_Strain_4,6),round(True_Plastic_Strain_5,6))

#print(round(Engg_Strain_1,6),round(Engg_Strain_2,6),round(Engg_Strain_3,6),round(Engg_Strain_4,6),round(Engg_Strain_5,6))

#print(Engg_Stress_1,Engg_Stress_2,Engg_Stress_3,Engg_Stress_4,Engg_Stress_5)

#print(True_Stress_1,True_Stress_2,True_Stress_3,True_Stress_4,True_Stress_5)

#Engg_Stress
Engg_Stress_6 = Engg_Stress_5 + 1*MyFty*0.05
Engg_Stress_7 = Engg_Stress_5 + 2*MyFty*0.05
Engg_Stress_8 = Engg_Stress_5 + 3*MyFty*0.05
Engg_Stress_9 = Engg_Stress_5 + 4*MyFty*0.05

#Engg_Strain
Engg_Strain_6 = Engg_Stress_6/myE+0.002*(Engg_Stress_6/MyFty)**Myn
Engg_Strain_7 = Engg_Stress_7/myE+0.002*(Engg_Stress_7/MyFty)**Myn
Engg_Strain_8 = Engg_Stress_8/myE+0.002*(Engg_Stress_8/MyFty)**Myn
Engg_Strain_9 = Engg_Stress_9/myE+0.002*(Engg_Stress_9/MyFty)**Myn

#True_Stress
True_Stress_6 = Engg_Stress_6*(1+Engg_Strain_6)
True_Stress_7 = Engg_Stress_7*(1+Engg_Strain_7)
True_Stress_8 = Engg_Stress_8*(1+Engg_Strain_8)
True_Stress_9 = Engg_Stress_9*(1+Engg_Strain_9)

#True_Plastic_Strain
True_Plastic_Strain_6 = np.log(1+Engg_Strain_6) - True_Stress_6/myE
True_Plastic_Strain_7 = np.log(1+Engg_Strain_7) - True_Stress_7/myE
True_Plastic_Strain_8 = np.log(1+Engg_Strain_8) - True_Stress_8/myE
True_Plastic_Strain_9 = np.log(1+Engg_Strain_9) - True_Stress_9/myE

#print(round(True_Plastic_Strain_6,6),round(True_Plastic_Strain_7,6),round(True_Plastic_Strain_8,6),round(True_Plastic_Strain_4,6),round(True_Plastic_Strain_9,6))
#print(True_Stress_8)

#print(round(Engg_Strain_6,6),round(Engg_Strain_7,6),round(Engg_Strain_8,6),round(Engg_Strain_4,6),round(Engg_Strain_9,6))

#print(Engg_Stress_6,Engg_Stress_7,Engg_Stress_8,Engg_Stress_9)
#print(True_Stress_6,True_Stress_7,True_Stress_8,True_Stress_9)

#Engg_Stress
Engg_Stress_10 = Engg_Stress_9 + 1*0.1*(MyFtu-MyFty)
Engg_Stress_11 = Engg_Stress_9 + 2*0.1*(MyFtu-MyFty)
Engg_Stress_12 = Engg_Stress_9 + 3*0.1*(MyFtu-MyFty)
Engg_Stress_13 = Engg_Stress_9 + 4*0.1*(MyFtu-MyFty)
Engg_Stress_14 = Engg_Stress_9 + 5*0.1*(MyFtu-MyFty)
Engg_Stress_15 = Engg_Stress_9 + 6*0.1*(MyFtu-MyFty)
Engg_Stress_16 = Engg_Stress_9 + 7*0.1*(MyFtu-MyFty)
Engg_Stress_17 = Engg_Stress_9 + 8*0.1*(MyFtu-MyFty)
Engg_Stress_18 = Engg_Stress_9 + 9*0.1*(MyFtu-MyFty)
Engg_Stress_19 = Engg_Stress_9 + 10*0.1*(MyFtu-MyFty)
Engg_Stress_20 = Engg_Stress_19

#Engg_Strain
Engg_Strain_10 = Engg_Stress_10/myE+0.002*(Engg_Stress_10/MyFty)**Myn
Engg_Strain_11 = Engg_Stress_11/myE+0.002*(Engg_Stress_11/MyFty)**Myn
Engg_Strain_12 = Engg_Stress_12/myE+0.002*(Engg_Stress_12/MyFty)**Myn
Engg_Strain_13 = Engg_Stress_13/myE+0.002*(Engg_Stress_13/MyFty)**Myn
Engg_Strain_14 = Engg_Stress_14/myE+0.002*(Engg_Stress_14/MyFty)**Myn
Engg_Strain_15 = Engg_Stress_15/myE+0.002*(Engg_Stress_15/MyFty)**Myn
Engg_Strain_16 = Engg_Stress_16/myE+0.002*(Engg_Stress_16/MyFty)**Myn
Engg_Strain_17 = Engg_Stress_17/myE+0.002*(Engg_Stress_17/MyFty)**Myn
Engg_Strain_18 = Engg_Stress_18/myE+0.002*(Engg_Stress_18/MyFty)**Myn
Engg_Strain_19 = Engg_Stress_19/myE+0.002*(Engg_Stress_19/MyFty)**Myn
Engg_Strain_20 = MyScu/100


#True_Stress
True_Stress_10 = Engg_Stress_10*(1+Engg_Strain_10)
True_Stress_11 = Engg_Stress_11*(1+Engg_Strain_11)
True_Stress_12 = Engg_Stress_12*(1+Engg_Strain_12)
True_Stress_13 = Engg_Stress_13*(1+Engg_Strain_13)
True_Stress_14 = Engg_Stress_14*(1+Engg_Strain_14)
True_Stress_15 = Engg_Stress_15*(1+Engg_Strain_15)
True_Stress_16 = Engg_Stress_16*(1+Engg_Strain_16)
True_Stress_17 = Engg_Stress_17*(1+Engg_Strain_17)
True_Stress_18 = Engg_Stress_18*(1+Engg_Strain_18)
True_Stress_19 = Engg_Stress_19*(1+Engg_Strain_19)
True_Stress_20 = Engg_Stress_20*(1+Engg_Strain_20)

#True_Plastic_Strain
True_Plastic_Strain_10 = np.log(1+Engg_Strain_10) - True_Stress_10/myE
True_Plastic_Strain_11 = np.log(1+Engg_Strain_11) - True_Stress_11/myE
True_Plastic_Strain_12= np.log(1+Engg_Strain_12) - True_Stress_12/myE
True_Plastic_Strain_13= np.log(1+Engg_Strain_13) - True_Stress_13/myE
True_Plastic_Strain_14= np.log(1+Engg_Strain_14) - True_Stress_14/myE
True_Plastic_Strain_15= np.log(1+Engg_Strain_15) - True_Stress_15/myE
True_Plastic_Strain_16= np.log(1+Engg_Strain_16) - True_Stress_16/myE
True_Plastic_Strain_17= np.log(1+Engg_Strain_17) - True_Stress_17/myE
True_Plastic_Strain_18= np.log(1+Engg_Strain_18) - True_Stress_18/myE
True_Plastic_Strain_19= np.log(1+Engg_Strain_19) - True_Stress_19/myE
True_Plastic_Strain_20= np.log(1+Engg_Strain_20) - True_Stress_20/myE


#print(round(True_Plastic_Strain_10,6),round(True_Plastic_Strain_11,6),round(True_Plastic_Strain_12,6),round(True_Plastic_Strain_13,6),round(True_Plastic_Strain_14,6),round(True_Plastic_Strain_15,6),round(True_Plastic_Strain_16,6),round(True_Plastic_Strain_17,6),round(True_Plastic_Strain_18,6),round(True_Plastic_Strain_19,6))
#print(round(Engg_Strain_10,6),round(Engg_Strain_11,6),round(Engg_Strain_12,6),round(Engg_Strain_13,6),round(Engg_Strain_14,6),round(Engg_Strain_15,6),round(Engg_Strain_16,6),round(Engg_Strain_17,6),round(Engg_Strain_18,6),round(Engg_Strain_19,6))


#print(Engg_Stress_10,Engg_Stress_11,Engg_Stress_12,Engg_Stress_13,Engg_Stress_14,Engg_Stress_15,Engg_Stress_16,Engg_Stress_17,Engg_Stress_18,Engg_Stress_19)
#print(True_Stress_10,True_Stress_11,True_Stress_12,True_Stress_13,True_Stress_14,True_Stress_15,True_Stress_16,True_Stress_17,True_Stress_18,True_Stress_19)
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

# variables

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

#Create ARC Point Solid line Of Plate

Create_ARC_Point_Solid_line_Create_Plate(myString,myPart_1,myArcX,myArcY,myArcDia,myEndLength,myPlateHalfHeight,myPlateLength,myPlateHalfThickness)

#------------------------------------------------------------------------------

#Create Bolt
Create_Bolt(myString,myPart_2,myArcX, myArcY,myBoltDia,myBoltHalfLength)

#------------------------------------------------------------------------------

#Plate material
Plate_material(myString,myMaterialName_1,myDensity,myE,myPoisonRatio)

#------------------------------------------------------------------------------

#Bolt materials
Bolt_materials(myString,myMaterialName_2,myDensity,myE,myPoisonRatio)

#------------------------------------------------------------------------------

#Create Section
Create_Section(myString,"Plate CS",myMaterialName_1)
Create_Section(myString,"Bolt CS",myMaterialName_2)

#------------------------------------------------------------------------------
#Create_Datum_Plane
myID_0 = Create_Datum_Plane(YZPLANE,myPart_1,myString,0.0)
myID_1 = Create_Datum_Plane(YZPLANE,myPart_2,myString,0.0)
myID_2 = Create_Datum_Plane(XZPLANE,myPart_2,myString,0.0)

#------------------------------------------------------------------------------

#Create_Partion
Create_Partion(myString,myPart_1,myID_0)
Create_Partion(myString,myPart_2,myID_1)
Create_Partion(myString,myPart_2,myID_2)

#------------------------------------------------------------------------------
#Section Assignment
Section_Assignment(myString,myPart_1,"Plate_1","Plate CS")
Section_Assignment(myString,myPart_2,"Bolt_1","Bolt CS")

#------------------------------------------------------------------------------

# Assemply
Assemply(myString,myPart_1,"Plate",0,0,0)
Assemply(myString,myPart_2,"Bolt",-myClearance/2,0, 0)

#------------------------------------------------------------------------------
# Reference point
myRP1,myRP_Position1 = Create_Reference_Point(-myClearance/2,0,0,myString,'RP-1')

#------------------------------------------------------------------------------

#RP to Bolt Rigid Body
Create_Interaction_Coupling(myString,"Bolt","RP-1","Bolt_1","RP_to_Bolt")

#------------------------------------------------------------------------------

#Create Step
Create_Step(myString,"Initial","Loading",ON,1500,0.01,1E-15,0.02)

#------------------------------------------------------------------------------

#FieldoutPut History Output Request
Fieldoutput_History_Output_Request(myString,"F-Output-1","H-Output-1",0.02)

#------------------------------------------------------------------------------

#Contact Property
Contact_Property(myString,"Intprop-1",0.3)

#------------------------------------------------------------------------------

#Bolt to Plate Arc Contact surface to surface

Bolt_Plate_Surface(myString,'Bolt','Plate','BS-1','PS-1')

#Bp_Contact(myString,"Bolt",'B_surface',"Plate",'A_surface','Int-1',"Loading","Intprop-1")

#------------------------------------------------------------------------------

#Create Contact Plate And Bolt Surface

Create_Contact_BP_Surface(myString,'BS-1','PS-1','Int-1','Loading','Intprop-1')

#------------------------------------------------------------------------------
# Bolt dsplacement boundary Conditions

Bolt_Displacement_Boundary_Conditions(myString,"RP-1",'Bolt_Displacement','Loading',-myDisplacement,0,0,0,0,0)

#------------------------------------------------------------------------------

Create_Surface(myString,myPart_1,myPlateLength-myEndLength,myPlateHalfHeight,myPlateHalfThickness,"Plate_Edge")
Create_Surface(myString,myPart_1,myPlateLength-myEndLength,0,myPlateHalfThickness/2,"Plate_Bottom_1")
Create_Surface(myString,myPart_1,-myEndLength,0,myPlateHalfThickness/2,"Plate_Bottom_2")
Create_Surface(myString,myPart_1,-myEndLength,myPlateHalfHeight/2,myPlateHalfThickness/2,"Plate_Edge_End_length")

#------------------------------------------------------------------------------

#Rigid Boundary Fixed Edged of plate
Create_BC_Plate_Edge_Fixed(myString,'Plate','Plate_Edge','Plate_Edge_Fixed','Loading')

#------------------------------------------------------------------------------

#Create_Extra_BC_For_Thin_Plate_Edge_Zaxis_Fixed(myString,'Plate','Plate_Edge_Fixed_End_L','Z_axis_fixed',myPlateLength/2,myPlateHalfHeight/2,myPlateHalfThickness)

#------------------------------------------------------------------------------

#Create Symmetry
Create_Symmetry(myString,'Plate','Plate_B_Symmetry','BC_Symmetry','Initial')

#------------------------------------------------------------------------------
Create_Z_Symmetry(myString,'Plate','Bolt','Z sym surf','Z Symmetry')
#------------------------------------------------------------------------------
#Create Mesh
Create_Mesh(myString,'Plate','Bolt',1.0,0.75)

#------------------------------------------------------------------------------

#Element Type
Select_Element_Type(myString,'Plate','Bolt',OFF,0.99)

#------------------------------------------------------------------------------

Create_Node_Set_ByBoundingBox(myString,'Plate',-myEndLength,0,0,0,myPlateHalfHeight,myPlateHalfThickness,"Element_NodeSet_1")

#------------------------------------------------------------------------------

Create_Element_Set_ByBoundingBox(myString,'Plate',-myEndLength,0,0,0,myPlateHalfHeight,myPlateHalfThickness,"Element_ElementSet_1")

#------------------------------------------------------------------------------

Create_Job(myString,myJobName,1)

#------------------------------------------------------------------------------

#os.chdir(myDirectory)

#------------------------------------------------------------------------------ 
#------------------------------------------------------------------------------ 

Create_Inp_File(myString)
#------------------------------------------------------------------------------

#Job submit
#Job_Submit(myJobName)

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
del mdb.models['Model-1']
#------------------------------------------------------------------------------

#Cae_File_Save(myfilesave) 
#mdb.jobs[myJobName].writeInput(consistencyChecking=OFF)

#------------------------------------------------------------------------------
