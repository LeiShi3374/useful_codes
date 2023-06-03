import os
import vtk
import numpy as np
import matplotlib.pyplot as plt

# User defined inputs
MESH_FILE = "../mesh/P1_0.005_Exact/mesh-complete.mesh.vtu"
RESULTS_FOLDER = "../results_mks_ma"
START_TIME = 1
END_TIME  = 100
DIFF_TIME = 1
TIME_STEP_SIZE = 0.01
PROBE_1 = [0.025, 0.03, 0.0]
PROBE_2 = [0.0, 0.03, 0.0]

#----------------------------------------------------------------------
def load_vtu(file_name):

    print "   Loading vtu file   <---   %s" % (file_name)
    vtuReader = vtk.vtkXMLUnstructuredGridReader()
    vtuReader.SetFileName(file_name)
    vtuReader.Update()

    vtuMesh = vtk.vtkUnstructuredGrid()
    vtuMesh = vtuReader.GetOutput()

    return vtuMesh
#----------------------------------------------------------------------

#----------------------------------------------------------------------
def get_closest_node(msh, xp):
    msh_npts = msh.GetNumberOfPoints()
    min_dist = 999999
    for ipt in range(0, msh_npts):
        msh_x = msh.GetPoint(ipt)
        dist  = np.sqrt((msh_x[0]-xp[0])**2 \
              + (msh_x[1]-xp[1])**2 + (msh_x[2]-xp[2])**2)
        if min_dist > dist:
            min_dist = dist
            ipt_near = ipt
            xpt_near = msh_x

    return ipt_near, xpt_near
#----------------------------------------------------------------------

#----------------------------------------------------------------------
def get_nodal_disps(file_name, ipt):
    msh = load_vtu(file_name)
    msh_disp = vtk.vtkDoubleArray()
    msh_disp = msh.GetPointData().GetArray("Displacement")

    return msh_disp.GetTuple3(ipt)
#----------------------------------------------------------------------

#----------------------------------------------------------------------
def write_disp_to_file(msh, pt, fname):

    (ipt, pcoords)  = get_closest_node(msh, pt)

    n_res = int((END_TIME - START_TIME)/DIFF_TIME) + 1
    u_res = np.zeros((n_res,3))
    time  = np.zeros((n_res,1))

    # Write displacements to file
    fid = open(fname,"w")
    fid.write("Variables = t, u_x, u_y, u_z\n")
    for i in range(0, n_res):
        ntime = START_TIME + i*DIFF_TIME
        time[i]  = float(ntime)*TIME_STEP_SIZE
        if (ntime < 100):
            vtu_file = "%s/result_%03d.vtu" %(RESULTS_FOLDER, ntime)
        else:
            vtu_file = "%s/result_%d.vtu" %(RESULTS_FOLDER, ntime)
        u_res[i,:] = get_nodal_disps(vtu_file, ipt)
        fid.write("%9.4f   %15.6e   %15.6e   %15.6e\n" % \
            (time[i], u_res[i,0], u_res[i,1], u_res[i,2]))
    fid.close()

    # Plot displacements using matplotlib
    plt.figure(figsize=(12,8))
    plt.plot(time, u_res[:,0], 'k-', label='u_x', linewidth=3.0)
    plt.plot(time, u_res[:,1], 'r--', label='u_y', linewidth=3.0)
    plt.plot(time, u_res[:,2], 'b-.', label='u_z', linewidth=3.0)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    axes = plt.gca()
    str_title = "Displacement vs. time @(%.3f %.3f %.3f)" % \
        (pt[0], pt[1], pt[2])
    axes.set_title(str_title)
    axes.set_xlabel("Time")
    axes.set_ylabel("Displacement (m)")

    axes.title.set_size(24)
    axes.xaxis.label.set_size(20)
    axes.yaxis.label.set_size(20)

    axes.grid(True)
    axes.legend(fontsize=18)

    fhdr, fext = os.path.splitext(fname)
    plt.savefig('fig_%s.png' % (fhdr))
    plt.show()


#----------------------------------------------------------------------

#----------------------------------------------------------------------
# Main function
if __name__ == '__main__':

    print "========================================================"

    # Load vtu mesh and find the nearest node
    msh = load_vtu(MESH_FILE)

    # Write displacement field for probe 1
    write_disp_to_file(msh, PROBE_1, "disp_ma_p1.dat")

    # Write displacement field for probe 1
    write_disp_to_file(msh, PROBE_2, "disp_ma_p2.dat")

    print "========================================================"

#----------------------------------------------------------------------

