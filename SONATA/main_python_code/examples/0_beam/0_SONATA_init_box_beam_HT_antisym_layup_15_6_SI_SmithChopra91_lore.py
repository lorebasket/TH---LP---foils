#  -*- coding: utf-8 -*-
"""
Created on Mo Oct 07 11:18:28 2019

@author: Roland Feil
"""

import os
import time
import numpy as np
import matplotlib.pyplot as plt

from SONATA.classBlade import Blade
from SONATA.utl.beam_struct_eval import beam_struct_eval


# ==== Main ==== #
start_time = time.time()

print('Current working directory is:', os.getcwd())


# ===== Provide Path Directory & Yaml Filename ===== #
folder_str = os.getcwd() + '/'
# run_dir = os.path.dirname( os.path.realpath(__file__) ) + os.sep
job_str = '0_box_beam_HT_antisym_layup_15_6_SI_SmithChopra91_lore.yaml'  # note: for better meshing convergence, units specified in yaml are in 'mm' instead of 'm'
job_name = 'box_beam_SmithChopra91_lore'
filename_str = folder_str + job_str


# ===== Define flags ===== #
flag_3d                 = False
flag_wt_ontology        = False
flag_ref_axes_wt        = False

# --- plotting fags --- #
mesh_resolution = 100              # it should be already defined in .yaml file
attribute_str           = 'MatID'   # default: 'MatID'
                                     # others:  'theta_3' - fiber orientation angle
                                     #          'stress.sigma11' (use sigma_ij to address specific component)                                            #          'stressM.sigma11'
                                     #          'strain.epsilon11' (use epsilon_ij to address specific component)                                           #          'strainM.epsilon11'
                                     #          'strainM.epsilon11'

# 2D cross sectional plats (blade_plot_sections)
flag_plotTheta11        = False     # plane orientation angle
flag_recovery           = False
flag_plotDisplacement   = True      # Needs recovery flag to be activated - shows displacements from loadings in cross sectional plots

# 3D plots (blade_post_3dtopo)
flag_wf                 = True      # plot wire-frame
flag_lft                = True      # plot lofted shape of blade surface (flag_wf=True obligatory); Note: create loft with grid refinement without too many radial_stations; can also export step file of lofted shape
flag_topo               = True      # plot mesh topology
c2_axis                 = False
flag_DeamDyn_def_transform            = False  # transform from SONATA to BeamDyn coordinate system
flag_write_BeamDyn = False                       # write BeamDyn input files for follow-up OpenFAST analysis (requires flag_DeamDyn_def_transform = True)
flag_write_BeamDyn_unit_convert = ''  #'mm_to_m'     # applied only when exported to BeamDyn files

# create flag dictionary
flags_dict = {"flag_wt_ontology": flag_wt_ontology, "flag_ref_axes_wt": flag_ref_axes_wt,
              "attribute_str": attribute_str,
              "flag_plotDisplacement": flag_plotDisplacement, "flag_plotTheta11": flag_plotTheta11,
              "flag_wf": flag_wf, "flag_lft": flag_lft, "flag_topo": flag_topo, "mesh_resolution": mesh_resolution,
              "flag_recovery": flag_recovery, "c2_axis": c2_axis}

# User defined radial stations
radial_stations = [0., 0.5, 1.]

# ==== Execute SONATA Blade Component Object ==== #
job = Blade(name=job_name, filename=filename_str, flags=flags_dict, stations=radial_stations)

# Build and mesh segments
job.blade_gen_section(topo_flag=True, mesh_flag = True, split_quads=True)

# ==== ANBAX Verification ==== #
# vabs_path = "/Users/rfeil/work/8_VABS/vabs_WIN/AnalySwift/VABS/VABSIII.exe"
# job.blade_run_vabs(vabs_path)

flag_3d                               = False
flag_run_vabs                         = False
flag_run_anbax                        = True
flag_verify_vabs_anbax                = False
flag_plot_vabs_struct_characteristics = False
flag_csv_export                       = True  # write BeamDyn input files for follow-up OpenFAST analysis (requires flag_DeamDyn_def_transform = True)

# Update flags dictionary
flags_dict ['flag_run_vabs'] = flag_run_vabs
flags_dict ['flag_run_anbax'] = flag_run_anbax
flags_dict ['flag_verify_vabs_anbax'] = flag_verify_vabs_anbax
flags_dict ['flag_DeamDyn_def_transform'] = flag_DeamDyn_def_transform
flags_dict ['flag_plot_vabs_struct_characteristics'] = flag_plot_vabs_struct_characteristics
flags_dict ['flag_csv_export'] = flag_csv_export
flags_dict ['flag_write_BeamDyn'] = flag_write_BeamDyn
flags_dict ['flag_write_BeamDyn_unit_convert'] = flag_write_BeamDyn_unit_convert
Loads_dict = {
	"Forces": [0., 0., 0.],  # [Fx, Fy, Fz]
	"Moments": [0., 0., 0.]   # [Mx, My, Mz]
}

beam_struct_eval(flags_dict=flags_dict, loads_dict=Loads_dict, vabs_path='', cs_pos=radial_stations, job=job, folder_str=folder_str, job_str=job_str, mu=None)

# ===== PLOTS ===== #
#job.blade_plot_attributes()
job.blade_plot_beam_props()

job.blade_plot_sections(attribute=attribute_str, plotTheta11=flag_plotTheta11, plotDisplacement=flag_plotDisplacement, savepath=folder_str)
if flag_3d:
    job.blade_post_3dtopo(flag_wf=flags_dict['flag_wf'], flag_lft=flags_dict['flag_lft'], flag_topo=flags_dict['flag_topo'])

print("--- Computational time: %s seconds ---" % (time.time() - start_time))

# EOF
