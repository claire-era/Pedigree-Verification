
This is a CLI based Java ARchive (.JAR) file for verifying progenies. This is made for the [International Rice Research Institute](http://irri.org/) for the internship project, <br><strong>Pedigree Verification.</strong>

[TOCM]

## Description
This project aims to identify which of the given F1 samples are true progenies.

## Dependencies
- This repository (Pedigree-Verification)
- `tassel-5-src`
                
1. Clone the repository by either calling `git clone https://github.com/InternationalRiceResearchInstitute/Pedigree-Verification` in the terminal <strong>OR</strong> clicking the `clone repository` button above.
2. Navigate to the location of the folder <strong>Pedigree-Verification</strong>.
##### Two Ways of Execution
                
    1. `java -jar <.hmp.txt> <pedigree_file.txt>`
    2. `java -jar <.hmp.txt> <pedigree_file.txt> <cut_off>`

    Or edit the script and run 
     `./start_default.sh` 
	 or 
	 `./start_cutoff.sh`

#####Notes
- `cut_off` should be of data type `double`.
- If the script requires permission to access, run `chmod 777 <script_name>`. This enables all file operations on the script for the current user.
                
##Limitations
#### File Formats
- Gene sequencing file - Only the HapMap text file is currently allowed in this package. There'll be future implementations of the project that can read other gene sequencing file formats such as: VCF, HDF5, etc.
- Pedigree file - The format of the pedigree file strictly should be as follows:

dnarun_name	 | germplasm_name | germplasm_pedigree | germplasm_type | germplasm_par1 | germplasm_par1_type |	germplasm_par2 | germplasm_par2_type | nasample_sample_group | dnasample_sample_group_cycle
------------- | -------------
- Only HapMap files(.hmp.txt) can be used. Again, future implementations will be done in the future.

##For Developers
The use of an IDE *(Recommended: Eclipse Photon to latest)* is highly recommended for developers to properly segregate the Java libraries.

<strong>Folder Structure</strong>
Note: When used in an IDE, put the repository inside the <strong>eclipse-workspace</strong> folder.
```
PedigreeVerification/
  bin/pedigreeVerification/
    MutableTaxaList.class
    PedVerification.class
    PedigreeFileInfo.class
  src/pedigreeVerification/
    MutableTaxaList.java
    PedVerification.java
    PedigreeFileInfo.java
  tassel5-src
  .settings/
  .classpath
  .project
  README.md
  qc2.geno.hmp.txt
  qc2.sample_no_spaces.txt
  PedVer.jar
  start_default.sh
  start_cutoff.sh
  
```
Where:
- `Pedver.jar` is the Runnable JAR to be executed by the scripts.
- `qc2.geno.hmp.txt` is a sample HapMap file
- `qc2.sample_no_spaces.txt` is a sample pedigree file
- `start_default.sh` is the default script for running default parameters 
- `start_cutoff.sh` is the script with indicated cut-off of scores

##About the Developer
This program is made by [Claire Era](https://github.com/claire-era), a student from University of the Philippines - Los Baños.

