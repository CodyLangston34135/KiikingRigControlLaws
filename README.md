# Kiiking Rig
Contains files required to run Kiiking rig (see example video) and process the data into angular displacement, velocity, and acceleration.

## Requirements
- MATLAB
- Image Processing Toolbox
- Signal Processing Toolbox

## Arduino Control Law
The main problem for the control laws is that the vibration from the motor adds lots of noise to the gyro at important times (peaks). 
The control law uses displacement/velocity data to find peaks in the data and attempts to backsolve the current frequency of the swing.
The control then attempts to move the sled up at the lowest point of the swing and down at the highest point of the swing.
