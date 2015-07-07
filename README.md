#![AFI icon](https://raw.githubusercontent.com/kjamison/afipeak/master/afipeak_icon_64x64.png) afipeak
MATLAB package for computing flip angles and suggested reference voltages from Dual-TR B1 mapping acquisitions, particularly for acquisitions with inhomogeneous B1 distribution.  

If you provide a maximum desired flip angle (default: 45&deg;+16%=52&deg;), this package will estimate the location of the highest flip angle in your AFI volume, and suggest a new reference voltage to produce the desired flip angle at that location.  The GUI provides a 3D volume browser to allow interactive placement of the "hot spot" ROI as well.

Flip angle computation from:
Yarnykh, V. L. (2007), Actual flip-angle imaging in the pulsed steady state: A method for rapid three-dimensional mapping of the transmitted radiofrequency field. Magn Reson Med, 57: 192–200. doi: 10.1002/mrm.21120

![AFI screenshot](https://raw.githubusercontent.com/kjamison/afipeak/master/afipeak_screenshot.png)

