import sys, rwkos

design_dir = rwkos.FindFullPath('siue/Research/SFLR_2010/vibration_suppression_design/TMM')

if design_dir not in sys.path:
    sys.path.append(design_dir)
