# MultifiberLocalization
MATLAB code to register CT scans to the Allen Mouse Brain CCF Atlas, and to localize multifiber photometry fibers. This code was written and developed in MATLAB2020b. See [Vu et al., 2023](https://www.biorxiv.org/content/10.1101/2023.11.17.567425v1).


Update 11/21/23: This repository is currently in progress and the code will be available shortly.



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
  * ../MRIAtlas/CCF/average_template_10_coronal.tif
  * ../MRIAtlas/CCF/annotation_10_coronal.tif
  * ../MRIAtlas/CCF/ccf_key.mat
  * ../MRIAtlas/CCF/landmarks.points
  * ../MRIAtlas/Chon/Chon_CCF_coronal.tif
  * ../MRIAtlas/Chon/Chon_labels_coronal.tif
  * ../MRIAtlas/Chon/chon_key.mat


![image](https://github.com/HoweLab/MultifiberLocalization/assets/21954946/3df9f697-58bc-45a8-b0d8-6db47faddd33)



## 2) REGISTER_CT_TO_ATLAS
Register a 3D CT image to the atlas.


## 3) LOCALIZE_FIBERS
Localize the fibers in a registered CT image.
   



  
# References   

## Code Dependencies 
The apps in this repository make use of the following repositories/tools. These have been packaged with the app. We note them here for reference.
* [point_to_line](https://github.com/thrynae/point_to_line_distance)
* [polyfitn](https://www.mathworks.com/matlabcentral/fileexchange/34765-polyfitn)
* [Fast_Tiff_Write](https://github.com/rharkes/Fast_Tiff_Write)
* [Bio-Formats](https://bio-formats.readthedocs.io/en/v7.0.1/users/matlab/index.html)
* [MATLAB Image Processing Toolbox](https://www.mathworks.com/products/image.html)


## Other
* https://github.com/cortex-lab/allenCCF: it is from here that we got the annotation labels for the atlases. 
  * We have created the following files:
    * CCF/ccf_key.mat  
    * Chon/chon_key.mat
  * Note that this repository also contains the very useful tool for Slice Histology Alignment, Registration, and Probe Track analysis: SHARP-Track.

