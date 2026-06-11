# IR View Suite

**WORK IN PROGRESS**

## Components

### IR View

TODO: Main application

### ircap

TODO: Server

### Other Utilities

The following global constants are defined in these utilities and control parameters such as resolution, frames per second and video device selection. Since camera characteristics vary between models, your mileage may vary:

- `WIDTH`: width of raw frame, in pixels
- `HEIGHT`: height of raw frame, in pixels
- `FRAME`: total frame size in bytes, derived from `WIDTH`, `HEIGHT`, and the number of bytes per pixel
- `MAX_FPS`: theoretical maximum frame rate of the camera
- `DEVICE`: path to the video device (usually `/dev/video*`)

#### irwebcam

Opens a camera stream via OpenCV viewer. Video is presented as received from the raw byte stream of the video device, with a color palette applied.

Since IR cameras typically provide relatively low-resolution images, the display can be enlarged using the `FACTOR` variable. When `ROI` is set to `True`, only a selected portion of the frame is displayed. The default region of interest targets the thermal image produced by the Hikmicro Mini2 camera, the device for which this utility was originally developed.

#### irshot

```sh
irshot [n_samples] [fps]
```

Captures a sequence of thermal frames from the camera and exports them to a MATLAB-compatible `.mat` file.

The first command-line argument specifies the number of frames to acquire, while the second specifies the desired sampling rate in frames per second. Since the camera operates at a fixed maximum frame rate, lower sampling rates are achieved by discarding intermediate frames. Thus, best consistency for delta time between frames is achieved by choosing a sampling rate that is easily divisible by your `MAX_FPS`.

Captured data are stored in a `.mat` file whose name includes the acquisition timestamp, number of samples, and sampling rate. The file contains the following arrays extracted from the raw camera stream:

- `TopImg`
- `MidLeftImg`
- `MidRightImg`
- `BottLeftImg`
- `BottRightImg`
- `Meta1`
- `Meta2`

These regions correspond to the different image and metadata areas present in the raw frame layout produced by the Hikvision Mini2 camera.

In addition to the image data, the output file stores a scalar variable named `frames_sec`, containing the sampling rate used during acquisition. This value is intended as input for IR View algorithms with an explicit dependence on the temporal spacing between samples (`Δt`).

## Dependencies

- `opencv-python`
- `scipy`
- `numpy`
- `matplotlib`
- `PyQt5`
