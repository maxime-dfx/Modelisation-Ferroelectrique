# -*- coding: utf-8 -*-

"""
Created on Tue Feb 24 08:23:23 2026

@author: lucas.lapostolle
"""




import numpy as np
from scipy.spatial import Voronoi


from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *


def csym(c11, c12, c13, c14, c15, c16, c22, c23, c24, c25, c26, c33, c34, c35, c36, c44, c45, c46, c55, c56, c66):
    return np.array([[c11, c12, c13, c14, c15, c16],
                     [c12, c22, c23, c24, c25, c26],
                     [c13, c23, c33, c34, c35, c36],
                     [c14, c24, c34, c44, c45, c46],
                     [c15, c25, c35, c45, c55, c56],
                     [c16, c26, c36, c46, c56, c66]])


def tens_4th_to_2nd(M):
    ans = np.zeros((6, 6))
    for i in range(6):
        for j in range(6):
            for k in range(3):
                for l in range(3):
                    for m in range(3):
                        for n in range(3):
                            ans[i, j] += beta(i+1, k+1, l+1) * \
                                M[k, l, m, n]*alpha(m+1, n+1, j+1)
    return ans


def tens_2nd_to_4th(M):
    ans = np.zeros((3, 3, 3, 3))
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    for m in range(6):
                        for n in range(6):
                            ans[i, j, k, l] += alpha(i+1, j+1,
                                                     m+1)*M[m, n]*beta(n+1, k+1, l+1)
    return ans


def alpha(i, j, k):
    ans = 0.
    if (i == j and j == k) and (i == 1 or i == 2 or i == 3):
        return 1.
    if (i == 2 and j == 3 and k == 4) or (i == 3 and j == 2 and k == 4) or (i == 1 and j == 3 and k == 5) or (i == 3 and j == 1 and k == 5) or (i == 1 and j == 2 and k == 6) or (i == 2 and j == 1 and k == 6):
        return 1/np.sqrt(2.)
    return ans

# @jit(nopython=True,cache=True)


def beta(i, j, k):
    ans = 0.
    if (i == j and j == k) and (i == 1 or i == 2 or i == 3):
        return 1.
    if (i == 4 and j == 2 and k == 3) or (i == 4 and j == 3 and k == 2) or (i == 5 and j == 1 and k == 3) or (i == 5 and j == 3 and k == 1) or (i == 6 and j == 1 and k == 2) or (i == 6 and j == 2 and k == 1):
        return 1/np.sqrt(2.)
    return ans


def rot_euler(phi1, phi, phi2):
    return np.array([[np.cos(phi1)*np.cos(phi2)-np.sin(phi1)*np.cos(phi)*np.sin(phi2), -np.cos(phi1)*np.sin(phi2)-np.sin(phi1)*np.cos(phi)*np.cos(phi2), np.sin(phi1)*np.sin(phi)],
                     [np.sin(phi1)*np.cos(phi2)+np.cos(phi1)*np.cos(phi)*np.sin(phi2), -np.sin(phi1)
                      * np.sin(phi2)+np.cos(phi1)*np.cos(phi)*np.cos(phi2), -np.cos(phi1)*np.sin(phi)],
                     [np.sin(phi)*np.sin(phi2), np.sin(phi)*np.cos(phi2), np.cos(phi)]])


def mat_rot_6(P):
    ans = np.zeros((6, 6), dtype=np.float32)
    for k in range(6):
        for l in range(6):
            for i in range(3):
                for j in range(3):
                    for m in range(3):
                        for n in range(3):
                            ans[k, l] += beta(k+1, i+1, j+1) * \
                                P[j, n]*P[i, m]*alpha(m+1, n+1, l+1)
    return ans


def rot_tensor_bis(C, P):
    P = mat_rot_6(P)
    ans = np.zeros((6, 6), dtype=np.float32)
    for i in range(6):
        for j in range(6):
            for k in range(6):
                for l in range(6):
                    ans[i, j] += P[i, k]*P[j, l]*C[k, l]
    return ans


def rotate_stiffness_tensor(c1111, c1122, c1212, angle):

    C = csym(c1111, c1122, c1122, 0, 0, 0, c1111, c1122, 0, 0, 0,
             c1111, 0, 0, 0, 2*c1212, 0, 0, 2*c1212, 0, 2*c1212)
    p = rot_euler(angle, 0, 0)
    C_rot = rot_tensor_bis(C, p)
    ans = tens_2nd_to_4th(C_rot)

    return (ans[0, 0, 0, 0], ans[0, 0, 1, 1], ans[1, 1, 1, 1], ans[0, 0, 2, 2], ans[1, 1, 2, 2], ans[2, 2, 2, 2], ans[0, 0, 0, 1], ans[1, 1, 0, 1], ans[2, 2, 0, 1], ans[0, 1, 0, 1], ans[0, 0, 0, 2], ans[1, 1, 0, 2], ans[2, 2, 0, 2], ans[0, 1, 0, 2], ans[0, 2, 0, 2], ans[0, 0, 1, 2], ans[1, 1, 1, 2], ans[2, 2, 1, 2], ans[0, 1, 1, 2], ans[0, 2, 1, 2], ans[1, 2, 1, 2])


if __name__ == "__main__":

    print("start")

    ############################################################
    # MODIFY THE FOLLOWING SECTION - Task 1
    ############################################################

    # Lx =     # dimension of the domain
    # Ly =     # dimension of the domain

    # c1111 =    # values of elastic constants 
    # c1122 = 
    # c1212 = 

    # N =    # number of grains

    # points =      # array containing the coordinates of seeds of voronoi tesselation, mirrored as in the preliminary section
                    # the creation of the voronoi tesselation is taken care automatically afterwards 

    # angles =      # array of random angles corresponding to the orientations of the grains

    ############################################################
    # DO NOT MODIFY AFTER THIS POINT
    ############################################################

    idx = np.argsort(points[:, 1])
    points = points[idx, :]

    voro = Voronoi(points, qhull_options='Qcc', furthest_site=False)
    vert = np.round(voro.vertices, 4)
    reg = voro.regions

    points = voro.points

    mdb.models['Model-1'].ConstrainedSketch(
        name='__profile__', sheetSize=200.0)
    mdb.models['Model-1'].sketches['__profile__'].rectangle(point1=(0.0-Lx*0.01, 0.0-Ly*0.01),
                                                            point2=(Lx+Lx*0.01, Ly+Ly*0.01))

    mdb.models['Model-1'].Part(dimensionality=TWO_D_PLANAR,
                               name='Part-1', type=DEFORMABLE_BODY)

    mdb.models['Model-1'].parts['Part-1'].BaseShell(
        sketch=mdb.models['Model-1'].sketches['__profile__'])

    mdb.models['Model-1'].parts['Part-1'].Set(
        faces=mdb.models['Model-1'].parts['Part-1'].faces.findAt(((Lx/2.0, Ly/2.0, 0),)), name='all')

    count = 0
    processed_vert = []

    facexleft = []
    facexright = []
    faceytop = []
    faceybottom = []

    for i in range(np.shape(points)[0]):
        if points[i, 0] >= 0 and points[i, 0] <= Lx and points[i, 1] >= 0 and points[i, 1] <= Ly:
            count += 1
            if count <= N-1:
                mdb.models['Model-1'].parts['Part-1'].DatumPointByCoordinate(
                    coords=(points[i, 0], points[i, 1], 0.0))
                grain = reg[voro.point_region[i]]

                for k in range(0, np.shape(grain)[0]):

                    if grain[k-1] != -1 and grain[k] != -1:

                        ############################################################
                        # MODIFY THE FOLLOWING SECTION - Task 2
                        ############################################################

                        # insert command line to create sktech line drawing the grains

                        ############################################################
                        # DO NOT MODIFY AFTER THIS POINT
                        ############################################################

                        if vert[grain[k-1], 0] == 0 and vert[grain[k], 0] == 0:
                            facexleft.append(
                                np.mean([vert[grain[k-1], 1], vert[grain[k], 1]]))
                        if vert[grain[k-1], 0] == Lx and vert[grain[k], 0] == Lx:
                            facexright.append(
                                np.mean([vert[grain[k-1], 1], vert[grain[k], 1]]))
                        if vert[grain[k-1], 1] == 0 and vert[grain[k], 1] == 0:
                            faceybottom.append(
                                np.mean([vert[grain[k-1], 0], vert[grain[k], 0]]))
                        if vert[grain[k-1], 1] == Ly and vert[grain[k], 1] == Ly:
                            faceytop.append(
                                np.mean([vert[grain[k-1], 0], vert[grain[k], 0]]))

                        processed_vert.append(
                            (vert[grain[k-1], 0], vert[grain[k-1], 1]))

                mdb.models['Model-1'].parts['Part-1'].PartitionFaceBySketch(faces=mdb.models['Model-1'].parts['Part-1'].faces.findAt(
                    ((points[i, 0], points[i, 1], 0),))[0], sketch=mdb.models['Model-1'].sketches['__profile__'])

            else:

                mdb.models['Model-1'].parts['Part-1'].DatumPointByCoordinate(
                    coords=(points[i, 0], points[i, 1], 0.0))
                grain = reg[voro.point_region[i]]

                for k in range(0, np.shape(grain)[0]):

                    if grain[k-1] != -1 and grain[k] != -1:

                        if vert[grain[k-1], 0] == 0 and vert[grain[k], 0] == 0:
                            facexleft.append(
                                np.mean([vert[grain[k-1], 1], vert[grain[k], 1]]))
                        if vert[grain[k-1], 0] == Lx and vert[grain[k], 0] == Lx:
                            facexright.append(
                                np.mean([vert[grain[k-1], 1], vert[grain[k], 1]]))
                        if vert[grain[k-1], 1] == 0 and vert[grain[k], 1] == 0:
                            faceybottom.append(
                                np.mean([vert[grain[k-1], 0], vert[grain[k], 0]]))
                        if vert[grain[k-1], 1] == Ly and vert[grain[k], 1] == Ly:
                            faceytop.append(
                                np.mean([vert[grain[k-1], 0], vert[grain[k], 0]]))
                            
            ############################################################
            # MODIFY THE FOLLOWING SECTION - Task 3
            ############################################################

            # insert command line to create the materials, the sections, and assign them to grains

            ############################################################
            # DO NOT MODIFY AFTER THIS POINT
            ############################################################

    del mdb.models['Model-1'].sketches['__profile__']
    mdb.models['Model-1'].ConstrainedSketch(gridSpacing=0.07, name='__profile__',
                                            sheetSize=2.88)
    mdb.models['Model-1'].parts['Part-1'].projectReferencesOntoSketch(
        filter=COPLANAR_EDGES, sketch=mdb.models['Model-1'].sketches['__profile__'])
    mdb.models['Model-1'].sketches['__profile__'].rectangle(point1=(0.0, 0.0),
                                                            point2=(Lx, Ly))
    mdb.models['Model-1'].sketches['__profile__'].rectangle(point1=(-Lx*0.01, -Ly*0.01),
                                                            point2=(Lx*(1+0.01), Ly*(1+0.01)))
    mdb.models['Model-1'].parts['Part-1'].Cut(
        sketch=mdb.models['Model-1'].sketches['__profile__'])

    mdb.models['Model-1'].rootAssembly.DatumCsysByDefault(CARTESIAN)
    mdb.models['Model-1'].parts['Part-1'].MaterialOrientation(
        additionalRotationType=ROTATION_NONE, axis=AXIS_3, fieldName='', localCsys=None, orientationType=GLOBAL, region=mdb.models['Model-1'].parts['Part-1'].sets['all'], stackDirection=STACK_3)
    mdb.models['Model-1'].rootAssembly.Instance(dependent=OFF, name='Part-1-1',
                                                part=mdb.models['Model-1'].parts['Part-1'])
    mdb.models['Model-1'].StaticStep(initialInc=0.01, maxInc=0.1, maxNumInc=10000,
                                     name='Step-1', previous='Initial')
    
    ############################################################
    # MODIFY THE FOLLOWING SECTION - Task 4
    ############################################################    
    
    
    # Create here the boundary conditions
    
    
    ############################################################
    # DO NOT MODIFY AFTER THIS POINT
    ############################################################
    

    mdb.models['Model-1'].rootAssembly.seedPartInstance(deviationFactor=0.1,
                                                        minSizeFactor=0.1, regions=(
                                                            mdb.models['Model-1'].rootAssembly.instances['Part-1-1'], ), size=max(Lx, Ly)/(N*100))
    mdb.models['Model-1'].rootAssembly.generateMesh(regions=(
        mdb.models['Model-1'].rootAssembly.instances['Part-1-1'], ))

    mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF,
            explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF,
            memory=90, memoryUnits=PERCENTAGE, model='Model-1', modelPrint=OFF,
            multiprocessingMode=DEFAULT, name='Voronoi_Job', nodalOutputPrecision=SINGLE,
            numCpus=1, numGPUs=0, numThreadsPerMpiProcess=1, queue=None, resultsFormat=ODB, scratch='', type=ANALYSIS, userSubroutine='', waitHours=0,
            waitMinutes=0)
    
    mdb.jobs['Voronoi_Job'].writeInput()

    ############################################################
    # MODIFY THE FOLLOWING SECTION - Task 5
    ############################################################