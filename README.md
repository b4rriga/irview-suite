# IR View Suite

**WORK IN PROGRESS**

## Components

### IR View

TODO: Main application

### ircap

TODO: Server

### Other Utilities

#### irwebcam

Open camera stream via OpenCV viewer. As IR cameras usually have a relatively low image resolution, a `FACTOR` variable may be modified, by which the resolution is multiplied. Video is presented as received from the raw byte stream of the video device, with a color palette applied.

![irwebcam](images/irwebcam.png)

#### irshot

## Dependencies

- `opencv-python`
- `scipy`
- `numpy`
- `matplotlib`
- `PyQt5`
