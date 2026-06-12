#!/usr/bin/env python3
#
# This file is part of IR View Suite.
#
# IR View Suite - Thermal/IR imaging suite
# Copyright (C) 2026 Hugo Barriga
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License v3
# as published by the Free Software Foundation.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY.
# See the GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU AGPLv3 along with this program.
# If not, see <https://www.gnu.org/licenses/>.

import configparser
import subprocess
import numpy as np
import cv2

cfg = configparser.ConfigParser()
cfg.read("/etc/irview.ini")

DEVICE = cfg["camera"]["device"]
WIDTH  = cfg.getint("camera", "width")
HEIGHT = cfg.getint("camera", "height")
FRAME  = WIDTH * HEIGHT * 2

FACTOR = 2
ROI    = False

p = subprocess.Popen(
	["v4l2-ctl", "-d", DEVICE, "--stream-mmap", "--stream-to=-", "--stream-count=0"],
	stdout=subprocess.PIPE,
	stderr=subprocess.DEVNULL
)

try:
	while True:
		frame = np.frombuffer(p.stdout.read(FRAME), dtype=np.uint16).reshape((HEIGHT, WIDTH))
		view = frame[196:292, 0:128] if ROI else frame
		norm = cv2.normalize(view, None, 0, 255, cv2.NORM_MINMAX, cv2.CV_8U)

		color = cv2.applyColorMap(norm, cv2.COLORMAP_INFERNO)

		size = (128 * FACTOR, 96 * FACTOR) if ROI else (WIDTH * FACTOR, HEIGHT * FACTOR)
		color = cv2.resize(color, (256*FACTOR, 344*FACTOR), interpolation=cv2.INTER_LINEAR)

		cv2.imshow("IR Webcam", color)
		if cv2.waitKey(1) == 27:
			break
finally:
	p.terminate()
	cv2.destroyAllWindows()
