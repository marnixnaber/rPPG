# Irissometry

Thank you for your interest in the rPPG toolbox for MATLAB!

This rPPG implementation does the following:
- Opens a popup window to select a directory with video files to be analyzed for heart rate
- Detects a person's face
- Detects and tracks unique face points to calculate face rotations over time
- Puts a mask on top of the face to detect only the skin (alternative color-based skin selection procedures are also available)
- Averages color changes across all facial skin pixels for a raw pulse signal
- Applies signal filtering and color space rotation algorithms to filter out movement and other noisy signals unrelated to heart rate pulsations (multiple algorithms available)
- Applies time frequency analysis to extract heart rate as a function of video time 
- Creates several plots to allow inspection of signal and HR detection results
- Saves the results per video file in an excel file

See "RunMe.m" for an example code. Just run it, select one or more videos, and observe the magic.

> Please cite our work in case you use our implementation!

![Example of face detection and tracking](https://github.com/marnixnaber/rPPG/blob/master/images/RPPG_image_output.png)

![Example of time frequency analysis](https://github.com/marnixnaber/rPPG/blob/master/images/RPPG_TFA.png)

![Example of heart rate as a function of video time](https://github.com/marnixnaber/rPPG/blob/master/images/RPPG_HR_time.png)

## Reference:

van der Kooij & Naber (2018). Standardized procedures for the testing and reporting of remote heart rate imaging. Behavior Research Methods

http://link.springer.com/article/10.3758/s13428-019-01256-8

If using the POS algorithm, please cite:
Plane-orthogonal-to-skin (pos) by Wang, W., den Brinker, A. C., Stuijk, S., & de Haan, G. (2017). Algorithmic principles of remote-PPG. IEEE. Transactions on Biomedical Engineering, 64(7), 1479-1491.

https://ieeexplore.ieee.org/abstract/document/7565547/

In case you plan to use rppg for  webcam videos collected online, then check out the following paper:
Di Lernia, D., Finotti, G., Tsakiris, M., Riva, G., & Naber, M. (2022, April 11). Remote photoplethysmography (rPPG) in the wild: Remote heart rate imaging via online webcams. https://doi.org/10.31234/osf.io/v89zn

https://psyarxiv.com/v89zn/

Also see the branch "into the wild" for an online data collection manual and adaptations to the rppg script to deal with settings and video database analyses more properly (thanks to Gianluca Finotti!).


## Info for videos:
Make sure that the videos have a high framerate (>20 frames per second) and lossless compression. Videos with VP8 and VP9 codecs are probably not supported and need to be converted to a lossless compression format (e.g., use ffmpeg). Resolution of the video is less relevant.


## Settings:

You can modify the parameters for the signal processing steps in the RunMe.m file.
For example, you can choose from multiple rPPG methods, including Naber (which involves an independent component analysis and several other processing steps) or the popular POS algorithm (see acknowledgments).

In the "extractFaceFromVideo.m" file you will find more parameters that can be set (e.g., sensitivity to detect faces, number of face points to track, and method to detect skin pixels from detected face).


## Required software:
The code has been tested in MATLAB 2019b on a Microsoft Windows 10 operating system. No guarantees can be provided for other MATLAB versions and operating platforms.

Please note that you probably need the following toolboxes to get the code to work:
- Computer vision toolbox 
- Image processing toolbox 
- Signal processing toolbox
- Statistics and machine learning toolbox

If you get an error about a missing function, it is likely that you have not installed the required toolboxes.

Please contact marnixnaber@gmail.com in case you cannot get the code to work. 
If you do so, please send along a screenshot of the error, the movie you want to analyze, and details regarding your operating system, matlab version, etc.

To get the scripts to work, download the following MATLAB scripts (also see acknowledgments):

- fastica (version 2.5): https://github.com/aludnam/MATLAB/tree/master/FastICA_25
- polyfitweighted: https://nl.mathworks.com/matlabcentral/fileexchange/13520-polyfitweighted



## License:

This updated work is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by-sa].

[![CC BY 4.0][cc-by-image]][cc-by-sa]

[cc-by-sa]: https://creativecommons.org/licenses/by-nc-sa/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by-nc-nd/4.0/88x31.png
<!-- https://i.creativecommons.org/l/by-nc-nd/4.0/88x31.png -->

Please contact marnixnaber@gmail.com for a commercial licence.


## Acknowledgments:

Plane-orthogonal-to-skin (pos) algorithm by Wang, W., den Brinker, A. C., Stuijk, S., & de Haan, G.

faceica matlab scripts by Hugo Gävert, Jarmo Hurri, Jaakko Särelä, and Aapo Hyvärinen.

polyfitweighted matlab script by Salman Rogers.

For more rPPG algorithms, see:

MATLAB code by Daniel McDuff's in github repository:
https://github.com/danmcduff/iphys-toolbox

PYTHON code by phuselab in github repository:
https://github.com/phuselab/pyVHR



## Contact:
For questions, please contact marnixnaber@gmail.com
