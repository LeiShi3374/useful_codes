import os
import sys
import vtk
import time
import numpy as np
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk

EPS = sys.float_info.epsilon
PI  = np.pi

DATASTR1 = 'Phi_EPI'
DATASTR2 = 'Phi_LV'
DATASTR3 = 'Phi_RV'
DATASTR4 = 'Phi_AB'

#----------------------------------------------------------------------
# Define thersholds for domains
# Epicardium - Myocardium
EPI_MYO = 0.8
# Endocardium - Myocardium
END_MYO = 0.8
# Medial - Basal
MID_BAS = 0.99
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# Define fiber and sheet angles
# Fiber angle at endocardium
ALFA_END = (+40.0)*PI/180.0
# Fiber angle at epicardium
ALFA_EPI = (-50.0)*PI/180.0
# Sheet angle at endocardium
BETA_END = (-65.0)*PI/180.0
# Sheet angle at epicardium
BETA_EPI = (+25.0)*PI/180.0
#----------------------------------------------------------------------

#----------------------------------------------------------------------
def isZero(a, b=None):

    absA = abs(a)
    if b == None:
        b = 0.0

    absA = abs(a)
    absB = abs(b)
    if absB > absA:
        absA, absB = absB, absA

    nrm  = max(absA, absB)
    flag = False
    if ((absA-absB)/nrm) < (10.0*EPS):
        flag = True

    return flag
#----------------------------------------------------------------------

#----------------------------------------------------------------------
def norm(u, v=None):

    if v == None:
        v = u

    n = np.size(u)
    l2norm = 0.0
    for i in range(0, n):
        l2norm = l2norm + (u[i]*v[i])

    return l2norm
#----------------------------------------------------------------------

#----------------------------------------------------------------------
def rot2quat(R):
    """
    ROT2QUAT - Transform Rotation matrix into normalized quaternion.
    Usage: q = rot2quat(R)
    Input:
    R - 3-by-3 Rotation matrix
    Output:
    q - 4-by-1 quaternion, with form [w x y z], where w is the scalar term.
    """
    # By taking certain sums and differences of the elements
    # of R we can obtain all products of pairs a_i a_j with
    # i not equal to j. We then get the squares a_i^2 from
    # the diagonal of R.
    a2_a3 = (R[0,1] + R[1,0]) / 4
    a1_a4 = (R[1,0] - R[0,1]) / 4
    a1_a3 = (R[0,2] - R[2,0]) / 4
    a2_a4 = (R[0,2] + R[2,0]) / 4
    a3_a4 = (R[1,2] + R[2,1]) / 4
    a1_a2 = (R[2,1] - R[1,2]) / 4
    D = np.array([[+1, +1, +1, +1],
                  [+1, +1, -1, -1],
                  [+1, -1, +1, -1],
                  [+1, -1, -1, +1]]) * 0.25
    aa = np.dot(D, np.r_[np.sqrt(np.sum(R**2) / 3), np.diag(R)])
    # form 4 x 4 outer product a \otimes a:
    a_a = np.array([[aa[0], a1_a2, a1_a3, a1_a4],
                 [a1_a2, aa[1], a2_a3, a2_a4],
                 [a1_a3, a2_a3, aa[2], a3_a4],
                 [a1_a4, a2_a4, a3_a4, aa[3]]])
    # use rank-1 approximation to recover a, up to sign.
    U, S, V = np.linalg.svd(a_a)
    q = U[:, 0]
    # q = np.dot(_math.sqrt(S[0]), U[:, 0])
    # Use this if you want unnormalized quaternions
    return q
#----------------------------------------------------------------------

#----------------------------------------------------------------------
def quat2rot(q):
    """
    QUAT2ROT - Transform quaternion into rotation matrix
    Usage: R = quat2rot(q)
    Input:
    q - 4-by-1 quaternion, with form [w x y z], where w is the scalar term.
    Output:
    R - 3-by-3 Rotation matrix
    """
    q = q / np.linalg.norm(q)
    w = q[0]; x = q[1];  y = q[2];  z = q[3]
    x2 = x*x;  y2 = y*y;  z2 = z*z;  w2 = w*w
    xy = 2*x*y;  xz = 2*x*z;  yz = 2*y*z
    wx = 2*w*x;  wy = 2*w*y;  wz = 2*w*z
    R = np.array([[w2+x2-y2-z2, xy-wz, xz+wy],
               [xy+wz, w2-x2+y2-z2, yz-wx],
               [xz-wy, yz+wx, w2-x2-y2+z2]])
    return R
#----------------------------------------------------------------------

#----------------------------------------------------------------------
def quat_mul(q0, q1):
    w0, x0, y0, z0 = q0
    w1, x1, y1, z1 = q1

    return np.array([-x1 * x0 - y1 * y0 - z1 * z0 + w1 * w0,
                     x1 * w0 + y1 * z0 - z1 * y0 + w1 * x0,
                     -x1 * z0 + y1 * w0 + z1 * x0 + w1 * y0,
                     x1 * y0 - y1 * x0 + z1 * w0 + w1 * z0])
#----------------------------------------------------------------------

#----------------------------------------------------------------------
def slerp(v0, v1, t):
    # >>> slerp([1,0,0,0],[0,0,0,1],np.arange(0,1,0.001))
    v0 = np.array(v0)
    v0 = v0/np.linalg.norm(v0)

    v1 = np.array(v1)
    v1 = v1/np.linalg.norm(v1)

    dot = np.sum(v0*v1)
    if (dot < 0.0):
        v1 = -v1
        dot = -dot

    DOT_THRESHOLD = 0.9995
    if (dot > DOT_THRESHOLD):
        result = v0+ t*(v1 - v0)
        result = result/np.linalg.norm(result)
        return result

    theta_0 = np.arccos(dot)
    sin_theta_0 = np.sin(theta_0)

    theta = theta_0*t
    sin_theta = np.sin(theta)

    s0 = np.cos(theta) - dot * sin_theta / sin_theta_0
    s1 = sin_theta / sin_theta_0

    return (s0 * v0)+ (s1 * v1)
#----------------------------------------------------------------------

#----------------------------------------------------------------------
def createVTKDataArray(dType, numComp, numTupls, dName):

    if dType == "double":
        D = vtk.vtkDoubleArray()
    elif dType == "int":
        D = vtk.vtkIntArray()

    D.SetNumberOfComponents(numComp)
    D.Allocate(numTupls)
    D.SetNumberOfTuples(numTupls)
    D.SetName(dName)

    return D
#----------------------------------------------------------------------

#----------------------------------------------------------------------
def loadLaplaceSoln(fileName):

    print "   Loading Laplace solution   <---   %s" % (fileName)
    vtuReader = vtk.vtkXMLUnstructuredGridReader()
    vtuReader.SetFileName(fileName)
    vtuReader.Update()

    vtuMesh = vtk.vtkUnstructuredGrid()

    print "   Extracting solution and its gradients at points"

    gradFilter = vtk.vtkGradientFilter()
    gradFilter.SetInputConnection(vtuReader.GetOutputPort())

    gradFilter.SetInputScalars(vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS,\
        DATASTR1)
    gradFilter.SetResultArrayName(DATASTR1 + '_pgrad')
    gradFilter.Update()
    vtuMesh = gradFilter.GetOutput()
    nPhiEP  = vtuMesh.GetPointData().GetArray(DATASTR1)
    nGPhiEP = vtuMesh.GetPointData().GetArray(DATASTR1+'_pgrad')

    gradFilter.SetInputScalars(vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS,\
        DATASTR2)
    gradFilter.SetResultArrayName(DATASTR2 + '_pgrad')
    gradFilter.Update()
    vtuMesh = gradFilter.GetOutput()
    nPhiLV  = vtuMesh.GetPointData().GetArray(DATASTR2)
    nGPhiLV = vtuMesh.GetPointData().GetArray(DATASTR2+'_pgrad')

    gradFilter.SetInputScalars(vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS,\
        DATASTR3)
    gradFilter.SetResultArrayName(DATASTR3 + '_pgrad')
    gradFilter.Update()
    vtuMesh = gradFilter.GetOutput()
    nPhiRV  = vtuMesh.GetPointData().GetArray(DATASTR3)
    nGPhiRV = vtuMesh.GetPointData().GetArray(DATASTR3+'_pgrad')

    gradFilter.SetInputScalars(vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS,\
        DATASTR4)
    gradFilter.SetResultArrayName(DATASTR4 + '_pgrad')
    gradFilter.Update()
    vtuMesh = gradFilter.GetOutput()
    nPhiAB  = vtuMesh.GetPointData().GetArray(DATASTR4)
    nGPhiAB = vtuMesh.GetPointData().GetArray(DATASTR4+'_pgrad')

    # Clean unnecessary arrays
    vtuMesh.GetPointData().RemoveArray(DATASTR1)
    vtuMesh.GetPointData().RemoveArray(DATASTR1+'_grad')

    vtuMesh.GetPointData().RemoveArray(DATASTR2)
    vtuMesh.GetPointData().RemoveArray(DATASTR2+'_grad')

    vtuMesh.GetPointData().RemoveArray(DATASTR3)
    vtuMesh.GetPointData().RemoveArray(DATASTR3+'_grad')

    vtuMesh.GetPointData().RemoveArray(DATASTR4)
    vtuMesh.GetPointData().RemoveArray(DATASTR4+'_grad')

    return vtuMesh, nPhiEP, nPhiLV, nPhiRV, nPhiAB, \
        nGPhiEP, nGPhiLV, nGPhiRV, nGPhiAB
#----------------------------------------------------------------------

#----------------------------------------------------------------------
def loadLaplaceSolnCells(fileName):

    print "   Loading Laplace solution at cells   <---   %s" % (fileName)
    vtuReader = vtk.vtkXMLUnstructuredGridReader()
    vtuReader.SetFileName(fileName)
    vtuReader.Update()

    vtuMesh = vtk.vtkUnstructuredGrid()

    print "   Extracting solution at cells"

    pt2Cell = vtk.vtkPointDataToCellData()
    pt2Cell.SetInputConnection(vtuReader.GetOutputPort())
    pt2Cell.PassPointDataOn()

    gradFilter = vtk.vtkGradientFilter()
    gradFilter.SetInputConnection(pt2Cell.GetOutputPort())

    gradFilter.SetInputScalars(vtk.vtkDataObject.FIELD_ASSOCIATION_CELLS,\
        DATASTR1)
    gradFilter.SetResultArrayName(DATASTR1 + '_grad')
    gradFilter.Update()
    vtuMesh = gradFilter.GetOutput()
    cPhiEP  = vtuMesh.GetCellData().GetArray(DATASTR1)
    cGPhiEP = vtuMesh.GetCellData().GetArray(DATASTR1+'_grad')

    gradFilter.SetInputScalars(vtk.vtkDataObject.FIELD_ASSOCIATION_CELLS,\
        DATASTR2)
    gradFilter.SetResultArrayName(DATASTR2 + '_grad')
    gradFilter.Update()
    vtuMesh = gradFilter.GetOutput()
    cPhiLV  = vtuMesh.GetCellData().GetArray(DATASTR2)
    cGPhiLV = vtuMesh.GetCellData().GetArray(DATASTR2+'_grad')

    gradFilter.SetInputScalars(vtk.vtkDataObject.FIELD_ASSOCIATION_CELLS,\
        DATASTR3)
    gradFilter.SetResultArrayName(DATASTR3 + '_grad')
    gradFilter.Update()
    vtuMesh = gradFilter.GetOutput()
    cPhiRV  = vtuMesh.GetCellData().GetArray(DATASTR3)
    cGPhiRV = vtuMesh.GetCellData().GetArray(DATASTR3+'_grad')

    gradFilter.SetInputScalars(vtk.vtkDataObject.FIELD_ASSOCIATION_CELLS,\
        DATASTR4)
    gradFilter.SetResultArrayName(DATASTR4 + '_grad')
    gradFilter.Update()
    vtuMesh = gradFilter.GetOutput()
    cPhiAB  = vtuMesh.GetCellData().GetArray(DATASTR4)
    cGPhiAB = vtuMesh.GetCellData().GetArray(DATASTR4+'_grad')

    # Clean unnecessary arrays
    vtuMesh.GetPointData().RemoveArray(DATASTR1)
    vtuMesh.GetCellData().RemoveArray(DATASTR1)
    vtuMesh.GetCellData().RemoveArray(DATASTR1+'_grad')

    vtuMesh.GetPointData().RemoveArray(DATASTR2)
    vtuMesh.GetCellData().RemoveArray(DATASTR2)
    vtuMesh.GetCellData().RemoveArray(DATASTR2+'_grad')

    vtuMesh.GetPointData().RemoveArray(DATASTR3)
    vtuMesh.GetCellData().RemoveArray(DATASTR3)
    vtuMesh.GetCellData().RemoveArray(DATASTR3+'_grad')

    vtuMesh.GetPointData().RemoveArray(DATASTR4)
    vtuMesh.GetCellData().RemoveArray(DATASTR4)
    vtuMesh.GetCellData().RemoveArray(DATASTR4+'_grad')

    return cPhiEP, cPhiLV, cPhiRV, cPhiAB
#----------------------------------------------------------------------

#----------------------------------------------------------------------
def setDomainID(vtuMesh, vtuPhiEP, vtuPhiLV, vtuPhiRV, vtuPhiAB):

    print "   Setting domain IDs at cells"
    numCell = vtuMesh.GetNumberOfCells()
    dmnIDs  = createVTKDataArray("int", 1, numCell, "DOMAIN_ID")

    for iCell in xrange(0, numCell):
        phiEP = vtuPhiEP.GetTuple1(iCell)
        phiLV = vtuPhiLV.GetTuple1(iCell)
        phiRV = vtuPhiRV.GetTuple1(iCell)
        phiAB = vtuPhiAB.GetTuple1(iCell)

        dmnID = 2
        if phiEP >= EPI_MYO:
            dmnID = 3
        if (phiLV >= END_MYO) | (phiRV >= END_MYO):
            dmnID = 1
        if phiAB >= MID_BAS:
            dmnID = 0
        dmnIDs.SetTuple1(iCell, dmnID)

    return dmnIDs
#----------------------------------------------------------------------

#----------------------------------------------------------------------
def getAlfaS(d):
    return ALFA_END*(1-d) - ALFA_END*d
#----------------------------------------------------------------------

#----------------------------------------------------------------------
def getBetaS(d):
    return BETA_END*(1-d) - BETA_END*d
#----------------------------------------------------------------------

#----------------------------------------------------------------------
def getAlfaW(d):
    return ALFA_END*(1-d) + ALFA_EPI*d
#----------------------------------------------------------------------

#----------------------------------------------------------------------
def getBetaW(d):
    return BETA_END*(1-d) + BETA_EPI*d
#----------------------------------------------------------------------

#----------------------------------------------------------------------
def axis (u, v):

    n1 = np.sqrt(norm(u))
    if not isZero(n1):
        e1 = u / n1
    else:
        e1 = [1.0, 0.0, 0.0]
        print " Zero gradient norm detected (e1)"

    e2 = v - np.dot(e1, v)*e1
    n2 = np.sqrt(norm(e2))
    if not isZero(n2):
        e2 = e2 / n2
    else:
        e2 = [0.0, 1.0, 0.0]
        print " Zero gradient norm detected (e2)"

    e0 = np.cross(e1, e2)

    Q  = np.zeros((3,3))
    Q[:,0] = e0
    Q[:,1] = e1
    Q[:,2] = e2

    return Q
#----------------------------------------------------------------------

#----------------------------------------------------------------------
def orient(Q, alpha, beta):

    ca = np.cos(alpha)
    sa = np.sin(alpha)
    cb = np.cos(beta)
    sb = np.sin(beta)

    Ra = np.array([ [ ca,  -sa,  0.0],
                    [ sa,   ca,  0.0],
                    [0.0,  0.0,  1.0]])

    Rb = np.array([ [1.0,  0.0,  0.0],
                    [0.0,   cb,   sb],
                    [0.0,  -sb,   cb]])

    Qt = np.matmul(Q, np.matmul(Ra, Rb))

    return Qt
#----------------------------------------------------------------------

#----------------------------------------------------------------------
def bislerp(Q_A, Q_B, t):

    qA = rot2quat(Q_A)
    qB = rot2quat(Q_B)

    i  = np.array([0, 1, 0, 0])
    j  = np.array([0, 0, 1, 0])
    k  = np.array([0, 0, 0, 1])

    qTest = np.zeros((4,8))
    qTest[:,0] =  qA
    qTest[:,1] = -qA
    qTest[:,2] =  quat_mul(i, qA)
    qTest[:,3] = -quat_mul(i, qA)
    qTest[:,4] =  quat_mul(j, qA)
    qTest[:,5] = -quat_mul(j, qA)
    qTest[:,6] =  quat_mul(k, qA)
    qTest[:,7] = -quat_mul(k, qA)

    norms=np.zeros(8)
    for i in xrange(0,8):
        norms[i] = np.linalg.norm(np.dot(qTest[:,i],qB))

    index = np.argmax(norms)
    qM = qTest[:,index]
    Q_AB = quat2rot(slerp(qM, qB, t))

    return Q_AB
#----------------------------------------------------------------------

#----------------------------------------------------------------------
def getFiberDirections(vtuMesh, vtuPhiEP, vtuPhiLV, vtuPhiRV, \
    vtuGPhiEP, vtuGPhiLV, vtuGPhiRV, vtuGPhiAB):

    numPoints = vtuMesh.GetNumberOfPoints()
    print "   Computing fiber directions at points"
    F = createVTKDataArray("double", 3, numPoints, "FIB_DIR1")
    S = createVTKDataArray("double", 3, numPoints, "FIB_DIR2")
    T = createVTKDataArray("double", 3, numPoints, "FIB_DIR3")

    j = 1
    k = 1
    print ("      Progress "),
    sys.stdout.flush()
    for iPt in xrange(0, numPoints):
        phiEP  = vtuPhiEP.GetTuple1(iPt)
        phiLV  = vtuPhiLV.GetTuple1(iPt)
        phiRV  = vtuPhiRV.GetTuple1(iPt)

        gPhiEP = vtuGPhiEP.GetTuple3(iPt)
        gPhiLV = vtuGPhiLV.GetTuple3(iPt)
        gPhiRV = vtuGPhiRV.GetTuple3(iPt)
        gPhiAB = vtuGPhiAB.GetTuple3(iPt)

        d = phiRV / max(EPS, phiLV + phiRV)
        alfaS = getAlfaS(d)
        betaS = getBetaS(d)
        alfaW = getAlfaW(phiEP)
        betaW = getBetaW(phiEP)

        Q_LV  = axis([gPhiAB[0], gPhiAB[1], gPhiAB[2]], \
            [-gPhiLV[0], -gPhiLV[1], -gPhiLV[2]])
        Q_LV  = orient(Q_LV, alfaS, betaS)

        Q_RV  = axis([gPhiAB[0], gPhiAB[1], gPhiAB[2]], \
            [ gPhiRV[0],  gPhiRV[1],  gPhiRV[2]])
        Q_RV  = orient(Q_RV, alfaS, betaS)

        Q_END = bislerp(Q_LV, Q_RV, d)

        Q_EPI = axis([gPhiAB[0], gPhiAB[1], gPhiAB[2]], \
            [ gPhiEP[0],  gPhiEP[1],  gPhiEP[2]])
        Q_EPI = orient(Q_EPI, alfaW, betaW)

        FST = bislerp(Q_END, Q_EPI, phiEP)

        F.SetTuple3(iPt, FST[0,0], FST[1,0], FST[2,0])
        S.SetTuple3(iPt, FST[0,1], FST[1,1], FST[2,1])
        T.SetTuple3(iPt, FST[0,2], FST[1,2], FST[2,2])
        if iPt==j:
            print ("%d%%  " % ((k-1)*10)),
            sys.stdout.flush()
            k = k + 1
            j = int(float((k-1)*numPoints)/10.0)
    print "[Done!]"

    return F, S, T
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# Main function
if __name__ == '__main__':

    t1 = time.time()
    print "========================================================"
    fileName = "laplace_soln.vtu"
    vtuMesh, phiEP, phiLV, phiRV, phiAB, \
        gPhiEP, gPhiLV, gPhiRV, gPhiAB = loadLaplaceSoln(fileName)

    F, S, T = getFiberDirections(vtuMesh, phiEP, phiLV, phiRV, \
        gPhiEP, gPhiLV, gPhiRV, gPhiAB)

    cphiEP, cphiLV, cphiRV, cphiAB = loadLaplaceSolnCells(fileName)
    dmnIDs = setDomainID(vtuMesh, cphiEP, cphiLV, cphiRV, cphiAB)

    print "   Writing fibers and domains to VTK data structure"
    vtuMesh.GetPointData().AddArray(F)
    vtuMesh.GetPointData().AddArray(S)
    vtuMesh.GetPointData().AddArray(T)

    vtuMesh.GetCellData().AddArray(dmnIDs)

    fileName = "fibers_domains.vtu"
    vtuWriter = vtk.vtkXMLUnstructuredGridWriter()
    vtuWriter.SetInputData(vtuMesh)
    vtuWriter.SetFileName(fileName)
    print "   Writing to vtu file   --->   %s" % (fileName)
    vtuWriter.Write()
    t2 = time.time()
    print('\n   Total time: %.3fs') % (t2-t1)
    print "========================================================"

#----------------------------------------------------------------------


