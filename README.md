# ISTTOK-GUI-PYTHON
GUI for ISTTok 
The files must be in the same folder.
Run the "main_gui_unina5_0.py" with python.
Sometimes the are problem with the "tomography", this isn't my function, to use it comment in "GUI_UNINA_5_0.py" the row 776-777 and uncomment 775.
  with tomography
  #tomography
            R_fromTomgraphy_value, z_fromTomgraphy_value = tomo_centroid(tomo_value)
            #R_fromTomgraphy_value = np.load("R_fromProbes.npy")
            #z_fromTomgraphy_value = np.load("z_fromProbes.npy")
            
  without tomogragphy
            #R_fromTomgraphy_value, z_fromTomgraphy_value = tomo_centroid(tomo_value)
            R_fromTomgraphy_value = np.load("R_fromProbes.npy")
            z_fromTomgraphy_value = np.load("z_fromProbes.npy"
