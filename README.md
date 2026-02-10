# MultifiberLocalization
MATLAB code to register CT scans to the Allen Mouse Brain CCF Atlas, and to localize multifiber photometry fibers. This code was written and developed in MATLAB2020b. See [Vu et al., 2024](https://www.sciencedirect.com/science/article/pii/S0896627323009704?via%3Dihub). The code used in the original publication has been improved and adapted for public use here. Contact maianhvu@bu.edu with any questions.



# Download
Before you begin, download the following program(s)/file(s).


## Fiji
https://imagej.net/software/fiji/




## Atlases: 
Download the appropriate files to generate the files needed for this pipeline:
1. Allen Mouse Brain CCF ([Wang et al., 2020](https://pubmed.ncbi.nlm.nih.gov/32386544/)):    
   * Download [average_template_10.nrrd](https://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/average_template/average_template_10.nrrd) 
   * Download [annotation_10.nrrd](https://download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/annotation/ccf_2022/annotation_10.nrrd)
   * More info can be found [here](https://help.brain-map.org/display/mouseconnectivity/API#API-DownloadAtlas3-DReferenceModels)
  
     
2. Kim Lab brain atlas ([Chon et al., 2019](https://pubmed.ncbi.nlm.nih.gov/31699990/)): 
   * Download [all data](https://kimlab.io/brain-map/atlas/assets/data_share/Atlas_Web_Release_data.7z) and unzip.
   * More info can be found [here](https://kimlab.io/brain-map/atlas/)    
  
# MAIN STEPS
Each of these steps is packaged in a MATLAB App. There are 2 ways to run each of these apps.
1. The .mlappinstall file contains the installer for the standalone app and includes the dependencies required.
2. If you have the dependencies required in your path (see below), you can alternatively just run the .mlapp file.


## 1) GENERATE_ATLAS_FILES 
This only needs to be run once to generate the atlas files necessary. Be sure you have atlas files (see above) downloaded first. You will end up with the following files in this repository folder (note that the double dots .. means one folder above this one):
  * **specified_atlas_location_path/MRIAtlas/CCF/average_template_10_coronal.tif** -- the atlas
  * **specified_atlas_location_path/MRIAtlas/CCF/annotation_10_coronal.tif** -- the annotated/labeled atlas
  * **specified_atlas_location_path/MRIAtlas/CCF/ccf_key.mat** -- key containing the annotation labels and other information
  * **specified_atlas_location_path/MRIAtlas/CCF/landmarks.points** -- the reference landmarks for registration
  * **specified_atlas_location_path/MRIAtlas/Chon/Chon_CCF_coronal.tif** -- the atlas
  * **specified_atlas_location_path/MRIAtlas/Chon/Chon_labels_coronal.tif** -- the annotated/labeled atlas
  * **specified_atlas_location_path/MRIAtlas/Chon/chon_key.mat** -- key containing the annotation labels and other information

![Screenshot 2023-12-19 145228](https://github.com/HoweLab/MultifiberLocalization/assets/21954946/1577bf6e-5e1f-4d1b-9298-072ae3d33b83)




## 2) REGISTER_CT
Register a 3D CT image to the atlas. See REGISTER_CT_README.pdf for instructions and more info.

![Screenshot 2023-12-19 151755](https://github.com/HoweLab/MultifiberLocalization/assets/21954946/1b18f51e-1816-43cf-a951-e9cd76d53725)




## 3) LOCALIZE_FIBERS
Localize the fibers in a registered CT image. See LOCALIZE_FIBERS_README.pdf for instructions and more info.

![Screenshot 2024-04-04 162019](https://github.com/HoweLab/MultifiberLocalization/assets/21954946/fe3d395c-eead-430e-bcf5-082a08814011)



## Other useful functions included:
1. **atlas_labels.m**
   * This function takes as input a m x 3 matrix of coordinates (mm from bregma) or indices (in atlas matrix space) in the order AP-ML-DV, where each row is a different point. It returns the anatomical labels assigned to the point by the Allen Mouse Brain CCF Atlas and the Kim Lab Atlas. Open the function or type 'help atlas_labels' for more info.
   * example: output = atlas_labels('coord', [0.8, 2.3, -3.4; 0.5, 2, -2]);
   * example: output = atlas_labels('idx', [445 800 413; 475 770 254]);
2. **howelab_table_labels.m**
   * Written for members of the Howe Lab, this function takes as input the path to the localized_fibers.mat struct, the struct itself, or the table contained in the "fiber_table" field of the struct, and returns a struct with a single field "table" which contains the same information as the fiber_table field but with fieldnames compatiable with Howe Lab legacy code. 
  
# References   

## Code Dependencies 
The apps in this repository make use of the following repositories/tools. These have been packaged with the app. We note them here for reference.
* [point_to_line](https://github.com/thrynae/point_to_line_distance)
* [polyfitn](https://www.mathworks.com/matlabcentral/fileexchange/34765-polyfitn)
* [Fast_Tiff_Write](https://github.com/rharkes/Fast_Tiff_Write)
* [Bio-Formats](https://bio-formats.readthedocs.io/en/v7.0.1/users/matlab/index.html)
* [MATLAB Image Processing Toolbox](https://www.mathworks.com/products/image.html)
* [MATLAB Statistics and Machine Learning Toolbox](https://www.mathworks.com/products/statistics.html)
* [sec2hms](https://www.mathworks.com/matlabcentral/fileexchange/22817-seconds-to-hours-minutes-seconds)


## Other
* https://github.com/cortex-lab/allenCCF: it is from here that we got the annotation labels for the atlases. 
  * We have created the following files:
    * CCF/ccf_key.mat  
    * Chon/chon_key.mat
  * Note that this repository also contains the very useful tool for Slice Histology Alignment, Registration, and Probe Track analysis: SHARP-Track.

------------------
##  License and Citation
 
If you use use this code in your research, please cite: 

Vu MT, Brown EH, Wen MJ, Noggle CA, Zhang Z, Monk KJ, Bouabid S, Mroz L, Graham BM, Zhuo Y, Li Y, Otchy TM, Tian L, Davison IG, Boas DA, Howe MW. [Targeted micro-fiber arrays for measuring and manipulating localized multi-scale neural dynamics over large, deep brain volumes during behavior. Neuron. 2024 Mar 20;112(6):909-923.e9.](https://www.sciencedirect.com/science/article/pii/S0896627323009704?via%3Dihub).

This repository is released under the [MIT License](https://opensource.org/license/mit) - see the [LICENSE](LICENSE) file for details.


## Acknowledgements

This work was supported by Aligning Science Across Parkinson’s **(ASAP-020370)** through the Michael J. Fox Foundation for Parkinson’s Research (MJFF), National Institute of Mental Health **(R01 MH125835)**, Whitehall Foundation Fellowship, Klingenstein-Simons Foundation Fellowship, and Parkinson’s Foundation (Stanley Fahn Junior Faculty Award, **PF-SF-JFA-836662**) to M.W.H.; NIMH **F32MH120894** to M.-A.T.V.
  


