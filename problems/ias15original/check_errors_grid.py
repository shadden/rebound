#!/usr/bin/env python

import subprocess
import os
import numpy as np
import matplotlib.pyplot as plt

ERRORS_FNAME = "errors_grid.txt"
N_MASS = 31
N_SIZE = 31
N_ORBITS = 10.
MASS_MIN = 1.e-3
MASS_MAX = 1.e3
SIZE_MIN = 1.e-3
SIZE_MAX = 1.e3

def call_rebound_grid(mass_scale_vec,size_scale_vec,ias15new=True):
    errors_A = np.NaN*np.zeros((len(mass_scale_vec),len(size_scale_vec)))
    # compile (make) binary
    version = "ias15"
    if ias15new:
        rebound_error_fname = "energy_ias15variable.txt"
    else:
        rebound_error_fname = "energy_ias15original.txt"
        version += "original"
    make_command = "make %s"%version
    subprocess.call(make_command,shell=True)
    errors_filename = ERRORS_FNAME[:-4] + '_' + version + ERRORS_FNAME[-4:]
    with open(errors_filename,'w') as errors_f:
        mass_scale_ndx = -1
        for mass_scale in mass_scale_vec:
            mass_scale_ndx += 1
            size_scale_ndx = -1
            for size_scale in size_scale_vec:
                size_scale_ndx += 1
                try: # remove the error file if it exists
                    os.remove(rebound_error_fname)
                except OSError: # swallow error and do nothing
                    pass
                runstring = "./nbody --scale_m=%8.4e --scale_a=%8.4e --tmax_factor=%4.2f"%(mass_scale,
                                                                                           size_scale,
                                                                                           N_ORBITS)
                print runstring
                subprocess.call(runstring,shell=True)
                # harvest errors from the rebound run
                rebound_errors = np.loadtxt(rebound_error_fname)
                this_error = rebound_errors[7]
                errors_f.write("%8.4e %8.4e %8.4e\n"%(mass_scale,
                                                      size_scale,
                                                      this_error))
                errors_A[mass_scale_ndx,size_scale_ndx] = this_error
    return errors_A

def generate_grid(ias15new=True):
    mass_scale_vec = np.logspace(np.log10(MASS_MIN),np.log10(MASS_MAX),N_MASS)
    size_scale_vec = np.logspace(np.log10(SIZE_MIN),np.log10(SIZE_MAX),N_SIZE)
    sizes_A,masses_A = np.meshgrid(size_scale_vec,mass_scale_vec)
    return dict(masses_A = masses_A,
                sizes_A  = sizes_A,
                errors_A = call_rebound_grid(mass_scale_vec,size_scale_vec,ias15new))

def read_grid(ias15new=True):
    version = "ias15"
    if ias15new: pass
    else: version += "original"
    errors_filename = ERRORS_FNAME[:-4] + '_' + version + ERRORS_FNAME[-4:]
    data = np.loadtxt(errors_filename)
    masses = data[:,0]
    sizes  = data[:,1]
    errors = data[:,2]
    new_shape = (N_MASS,N_SIZE)
    return dict(masses_A = np.reshape(masses,new_shape),
                sizes_A  = np.reshape(sizes,new_shape),
                errors_A = np.reshape(errors,new_shape))

def display_grid(info_dict,ias15new=True):
    masses_A = info_dict["masses_A"]
    sizes_A  = info_dict["sizes_A"]
    errors_A = info_dict["errors_A"]
    log_abs_errors_A = np.log10(np.abs(errors_A + 1.00001e-17))
    errors_A = log_abs_errors_A
    mass_vec = masses_A[:,0]
    size_vec = sizes_A[0,:]
    logmass_vec_edges = np.zeros(len(mass_vec)+1)
    logsize_vec_edges = np.zeros(len(mass_vec)+1)
    logmass_vec_edges[:-1] = np.log10(mass_vec) - 0.5*(np.log10(mass_vec[1])
                                                       - np.log10(mass_vec[0]))
    logmass_vec_edges[-1] = np.log10(mass_vec[-1]) + 0.5*(np.log10(mass_vec[1])
                                                          - np.log10(mass_vec[0]))
    mass_vec_edges = 10.**logmass_vec_edges
    logsize_vec_edges[:-1] = np.log10(size_vec) - 0.5*(np.log10(size_vec[1])
                                                       - np.log10(size_vec[0]))
    logsize_vec_edges[-1] = np.log10(size_vec[-1]) + 0.5*(np.log10(size_vec[1])
                                                          - np.log10(size_vec[0]))
    size_vec_edges = 10.**logsize_vec_edges
    plt.ion()
    plt.pcolor(size_vec_edges,mass_vec_edges,errors_A)
    plt.xlabel("Size Scale")
    plt.ylabel("Mass Scale")
    plt.gca().set_xscale('log')
    plt.gca().set_yscale('log')
    plt.colorbar()
    plt.autoscale(tight=True)
    if ias15new:
        plt.title("log10[|Errors|] (new time-stepping)")
        plt.show()
        plt.savefig("errors_new.pdf")
        plt.close()
    else:
        plt.title("log10[|Errors|] (original time-stepping)")
        plt.show()
        plt.savefig("errors_original.pdf")
        plt.close()
    return 1

if __name__=="__main__":
    ias15new_list = [False,True]
    #ias15new_list = [False]
    for ias15new_bool in ias15new_list:
        make_grid = True
        #make_grid = False
        if make_grid:
            ret_dict = generate_grid(ias15new=ias15new_bool)
        else: # read grid from file
            ret_dict = read_grid(ias15new=ias15new_bool)
        display_grid(ret_dict,ias15new=ias15new_bool)
