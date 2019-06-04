# rPPG - remote heart rate detection / photoplethysmography software development

run "RunMe.m" for a demo. Adapt this script to you own needs.

Run in matlab version 2014b or younger. 

To get the scripts to work, first download the fastica (version 2.5) matlab scripts
by Hugo Gävert, Jarmo Hurri, Jaakko Särelä, and Aapo Hyvärinen:

https://github.com/aludnam/MATLAB/tree/master/FastICA_25



% -----------------ACKNOWLEDGMENT-------------------
% van der Kooij & Naber (2018). Standardized procedures for the testing and 
% reporting of remote heart rate imaging. Behavior Research Methods
%
%
% -----------------CONTACT-------------------
% 
% For questions, remarks, or problems with the code, 
% please contact: marnixnaber@gmail.com
%
%
% -----------------PARAMETERS-------------------
%
% Below you can vary the parameters for the signal processing steps (e.g.
% frequency filtering).
% 
% In the "extractFaceFromVideo.m" file you will find more parameters that
% can be adjusted (e.g., sensitivity to detect faces, number of points to
% track the face, and method to detect skin pixels)
%

% --------LICENSE & ACKNOWLEDGMENT-----------
% 
% Copyright © 2017 Marnix Naber, The Netherlands
% 
% This program is distributed under the terms of the GNU General Public
% % License v3.0 (see gpl.txt)
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
% 
