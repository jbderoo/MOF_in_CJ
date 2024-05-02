![Logo](images/Graphical_Abstract.png "Graphical Abstract")
### Cu_in_pore_calcs: ###
  Computational determination of the maximum amount of CuBTC that can fit inside of a CJ nanopore. This is used to help compute the packing fraction as the upper limit number that could fit inside the pore.
  
  A_center_mof.py is used to take the tiled CuBTC MOF (in a 7x7x5 unit cell size) and center it. These files are cubtc_775_no_dups.pdb and cubtc_775_no_dups_centered.pdb, respectively.
  
  B_chainAllA.py is used to create a pdb file where all of the CJ monomers have the same chain. This is done to simplify the computation in later steps, where the pore and the MOF are used at the same time to detect clashes and remove them. Two files (221_fixed_CJ_centered_rechained_part1.pdb and 221_fixed_CJ_centered_rechained_part2.pdb) should be used for viewing. There are two files for viewing because there are more than 26 monomers present required to create the pore, but PyMOL cannot render two chain As independently, easily. Therefore part1 contains chains A-Z, and part2 contains chains A-J, to fully construct the pore. 221_CJ_centered.pdb is the combination of these two files via `cat 221_fixed_CJ_centered_rechained_part1.pdb 221_fixed_CJ_centered_rechained_part2.pdb > 221_CJ_centered.pdb`. The script B_chainAllA.py then ensures 221_CJ_centered.pdb only contains chain A.
  
  C_punchout.py is used to drastically cut down on the size of the MOF before doing clash detection calculations, as these are fairly computationally expensive. A radius is created where everything much larger than the pore is immediately thrown out.
  
  D_clash_detection_pore_cubtc.py then does clash detection and further refines random little edge effects of the MOF@CJ interaction. It also reports back the theoretical upper limit of Cu inside the nanopore, 2059.

### scXRD_to_PXRD: ###
  Example of how the data collected (images in results/) is turned into a pXRD plot for comparison. This specific example is the 90 degree rotation analysis. Included are images and scripts for doing the CuBTC and UiO-67 analysis in their respective directories. These scripts rely on the Fabio package (https://pypi.org/project/fabio/). I only use it to extract the data from the .img files, I could not get their automatic integration tools to work properly for me (likely user error). I wrote my own, they're found in xrd_tools.py, which is found in each of these subdirectories.

  To help explain the idea of azimuthal integration (especially if the paper description and still image wasn't clear), I've prepared this gif. As we start in the center of the image and walk outwards, we average the radius. We can the convert this distance from the center to a literal, physical distance, then use Brag's law and geometry to compute a 2theta. We then can compare this computed, extracted pXRD to the expected pXRD to verify our material of interest is actually present (CuBTC!). If the gif is stuck at the end, try reloading the page - sometimes it loops endlessly for me, and othertimes it doesn't.

![Azimuthal Integration](images/aziint_compressed.gif "Azimuthal Integration GIF")

