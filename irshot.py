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
import numpy as np

from datetime import datetime
from scipy.io import savemat

DEVICE = "/dev/video2"

if len(sys.argv) < 2:
    nframes = 1
else:
    if sys.argv[1] == "help" or sys.arvg[1] == "?":
        print(f"usage: {sys.argv[0]} [n_samples]") # TODO: add delta_t argument
    try:
        nframes = int(sys.argv[1])
        if nframes < 1:
            raise ValueError
    except ValueError:
        print("error: invalid number for amount of samples")
        sys.exit(1)

p = subprocess.Popen(
    ["v4l2-ctl", "-d", DEVICE, "--stream-mmap", "--stream-to=-", "--stream-count=0"],
    stdout=subprocess.PIPE,
    stderr=subprocess.DEVNULL
)

WIDTH  = 256
HEIGHT = 344
FRAME  = WIDTH * HEIGHT * 2

caps = {
    "TopImg":       np.empty((192, 256, nframes), dtype=np.uint16),
    "MidLeftImg":   np.empty((96,  128, nframes), dtype=np.uint16),
    "MidRightImg":  np.empty((96,  128, nframes), dtype=np.uint16),
    "BottLeftImg":  np.empty((48,  128, nframes), dtype=np.uint16),
    "BottRightImg": np.empty((48,  128, nframes), dtype=np.uint16),
    "Meta1":        np.empty((4,   256, nframes), dtype=np.uint16),
    "Meta2":        np.empty((4,   256, nframes), dtype=np.uint16)
}

try:
    for t in range(nframes):

        raw = np.frombuffer(
            p.stdout.read(FRAME),
            dtype=np.uint16
        )

        if raw.size != WIDTH * HEIGHT:
            raise RuntimeError(
                f"Uncomplete frame: expected "
                f"{WIDTH * HEIGHT} pixels, "
                f"received {raw.size}"
            )

        frame = raw.reshape((HEIGHT, WIDTH))

        caps["TopImg"][:, :, t]       = frame[0:192,   0:256]
        caps["MidLeftImg"][:, :, t]   = frame[196:292, 0:128]
        caps["MidRightImg"][:, :, t]  = frame[196:292, 128:256]
        caps["BottLeftImg"][:, :, t]  = frame[292:340, 0:128]
        caps["BottRightImg"][:, :, t] = frame[292:340, 128:256]
        caps["Meta1"][:, :, t]        = frame[192:196, 0:256]
        caps["Meta2"][:, :, t]        = frame[340:344, 0:256]
finally:
    p.terminate()
    p.wait()

stamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
savemat(stamp + f"_[{nframes}]" + "_irshot.mat", caps) # TODO: include delta_t into output file name
