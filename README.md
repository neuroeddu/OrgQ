# OrgQ

A Fiji macro written in jython to analysis IF images of cyrosection from organoid tissues. 

Rhalena Thomas, Eddie Cai and Harris Nami

We provide analysis 4 (or more) channels IF of histology sections of organoids.  Written for midbrain organoid images.

This is a macro written in Jython for Fiji (ImageJ) you must use Fiji because some 
Pull program into Fiji and press run.  To see the images and optimize settings set 'batch mode = FALSE'

You will be asked to select your input folder. This will contain images of organoids to quantify. Select your output folder.  
The results will be a csv file with data summarized by image.

The macro requires confocal images with file_name_CH00, CH01, etc as well as a merge image.  The prefix must match between each file name and spaces are not accepted. 
