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

import os
import sys
import subprocess
import configparser
import numpy as np

from datetime import datetime
from scipy.io import savemat

cfg = configparser.ConfigParser()
cfg.read("/etc/irview.ini")

DEVICE  = cfg["camera"]["device"]
WIDTH   = cfg.getint("camera", "width")
HEIGHT  = cfg.getint("camera", "height")
FRAME   = WIDTH * HEIGHT * 2
MAX_FPS = cfg.getint("camera", "max_fps")

def parse_args():
	if len(sys.argv) < 2:
		return 1, MAX_FPS

	if sys.argv[1] in ("help", "h", "?"):
		print(f"usage: {os.path.basename(sys.argv[0])} [n_samples] [fps]")
		sys.exit(0)

	try:
		nframes = int(sys.argv[1])
		if nframes < 1:
			raise ValueError
	except ValueError:
		print("error: invalid number of samples")
		sys.exit(1)

	fps = MAX_FPS
	if len(sys.argv) >= 3:
		try:
			fps = float(sys.argv[2])
			if fps <= 0 or fps > MAX_FPS:
				raise ValueError
		except ValueError:
			print("error: fps must be in (0, 25]")
			sys.exit(1)

	return nframes, fps

nframes, fps = parse_args()

decim = int(round(MAX_FPS / fps))
if decim < 1:
	decim = 1

p = subprocess.Popen(
	["v4l2-ctl", "-d", DEVICE, "--stream-mmap", "--stream-to=-", "--stream-count=0"],
	stdout=subprocess.PIPE,
	stderr=subprocess.DEVNULL
)

caps = {
	"frames_sec":   np.array(fps, dtype=np.float32),
	"TopImg":       np.empty((192, 256, nframes), dtype=np.uint16),
	"MidLeftImg":   np.empty((96,  128, nframes), dtype=np.uint16),
	"MidRightImg":  np.empty((96,  128, nframes), dtype=np.uint16),
	"BottLeftImg":  np.empty((48,  128, nframes), dtype=np.uint16),
	"BottRightImg": np.empty((48,  128, nframes), dtype=np.uint16),
	"Meta1":        np.empty((4,   256, nframes), dtype=np.uint16),
	"Meta2":        np.empty((4,   256, nframes), dtype=np.uint16)
}

try:
	t = 0
	grabbed = 0

	while grabbed < nframes:
		raw = np.frombuffer(p.stdout.read(FRAME), dtype=np.uint16)

		if raw.size != WIDTH * HEIGHT:
			continue

		if t % decim != 0:
			t += 1
			continue

		frame = raw.reshape((HEIGHT, WIDTH))

		caps["TopImg"][:, :, grabbed]       = frame[0:192,   0:256]
		caps["MidLeftImg"][:, :, grabbed]   = frame[196:292, 0:128]
		caps["MidRightImg"][:, :, grabbed]  = frame[196:292, 128:256]
		caps["BottLeftImg"][:, :, grabbed]  = frame[292:340, 0:128]
		caps["BottRightImg"][:, :, grabbed] = frame[292:340, 128:256]
		caps["Meta1"][:, :, grabbed]        = frame[192:196, 0:256]
		caps["Meta2"][:, :, grabbed]        = frame[340:344, 0:256]

		grabbed += 1
		t += 1

finally:
	p.terminate()
	p.wait()

stamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
out_name = f"{stamp}_[{nframes}]_[{fps:.2f}fps]_irshot.mat"

savemat(out_name, caps)
print(f"Saved `{os.path.abspath(out_name)}` with {nframes} frames captured at {fps} fps")
