**Fully automated pipeline for detection of sex linked genes using RNA-Seq data**

This pipeline has been designed for detection of sex-linked genes from the cross using RNASeq data. 
It contains two workflows, one for detection of Y-linked and one for detection of X-linked genes.

INSTALLING PREREQUISITES
You will need 3 prerequisites in order to successfully run this pipeline:
1) install linkyx
2) install tool dependencies
3) import workflows

1) The easiest way how to install linkyx tools is to login to your local Galaxy instance and click 
Admin -> Search and browse tool sheds -> Galaxy main tool shed -> linkyx

Click on the tool and then simply press button "Install to Galaxy" in the upper part of the screen.
This tool repository contains four tools:
LINKYX_X
Identify X linked contigs from cross

LINKYX_Y
Identify Y linked contigs from cross

LINKYX_mpileup
Run mpileup with LINKYX-specific parameters

Reformat Trinity header

All tools can be run independently, but the easiest way to obtain proper input files is via provided workflows.

2) The dependencies can be installed in a similar fashion as linkyx tool from the Galaxy tool shed by searching 
for suite_linkyx_1_0

3) When installing linkyx from toolshed, you can click on the "Workflows" on the main page describing the tool 
in the section "Contents of this repository". After clicking on either linkY or linkX, you will be presented an 
option to install the workflows by clicking "Install to Galaxy" button.

RUNNING THE WORKFLOWS
After you imported workflows to the galaxy and installed linkyx and workflow dependencies, the pipeline is ready to go. 

You can access it by clicking on the "Workflow" in the upper gray layer of your Galaxy instance. After clicking on the 
name of the workflow and choosing "Run", you will be requested to select input datasets, two for each member of the cross 
since paired-end reads are expected by default. If you don't have your data yet, test files are available in this repository. 
The computation for these files ends in ~5 minutes and outputs one X-linked and one Y-linked contig.

CONTACT INFORMATION
If you have any questions, please use following email address: biomonika@psu.edu