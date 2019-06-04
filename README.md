# rPPG - remote heart rate detection / photoplethysmography software development

run "RunMe.m" for a demo. Adapt this script to you own needs.

Run in matlab version 2014b or younger. 

To get the scripts to work, first download the fastica (version 2.5) matlab scripts
by Hugo Gävert, Jarmo Hurri, Jaakko Särelä, and Aapo Hyvärinen:

https://github.com/aludnam/MATLAB/tree/master/FastICA_25


-----------------ACKNOWLEDGMENT-------------------

van der Kooij & Naber (2018). Standardized procedures for the testing and reporting of remote heart rate imaging. Behavior Research Methods

http://link.springer.com/article/10.3758/s13428-019-01256-8


-----------------CONTACT-------------------
 
For questions, remarks, or problems with the code, please contact: marnixnaber@gmail.com


-----------------PARAMETERS-------------------

Below you can vary the parameters for the signal processing steps (e.g. frequency filtering).

In the "extractFaceFromVideo.m" file you will find more parameters that can be adjusted (e.g., sensitivity to detect faces, number of points to track the face, and method to detect skin pixels)
