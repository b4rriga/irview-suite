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

import socket
import struct
import time
import subprocess
import numpy as np
import threading
import os

_CAM_MAGIC   = b"IRFR"
_CAM_HDR_FMT = "<4sIIId"

DEVICE = "/dev/video2"

WIDTH  = 256
HEIGHT = 344
FRAME  = WIDTH * HEIGHT * 2

X0, X1 = 0, 256
Y0, Y1 = 0, 192
W = X1 - X0
H = Y1 - Y0

def tv_bars():
	f = np.zeros((H, W), dtype=np.float32)
	n = 8
	step = max(1, W // n)
	vals = np.linspace(0.1, 1.0, n)

	for i in range(n):
		f[:, i * step:(i + 1) * step] = vals[i]

	f[::2]  *= 0.6
	f[0, :]  = 1
	f[-1, :] = 1
	f[:, 0]  = 1
	f[:, -1] = 1
	return f

TV_BARS = tv_bars()

class Camera:
	def __init__(self):
		self.lock = threading.Lock()
		self.frame = None
		self.proc = None
		self.buf = bytearray()

		threading.Thread(target=self._loop, daemon=True).start()

	def _start(self):
		if not os.path.exists(DEVICE):
			return None

		return subprocess.Popen(
			["v4l2-ctl", "-d", DEVICE, "--stream-mmap", "--stream-to=-", "--stream-count=0"],
			stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, bufsize=0
		)

	def _loop(self):
		while True:
			if self.proc is None:
				self.proc = self._start()
				time.sleep(0.2)
				continue

			chunk = self.proc.stdout.read(4096)

			if not chunk:
				self.proc = None
				continue

			self.buf.extend(chunk)

			while len(self.buf) >= FRAME:
				raw = self.buf[:FRAME]
				del self.buf[:FRAME]

				frame = np.frombuffer(raw, dtype=np.uint16).reshape((HEIGHT, WIDTH))
				roi = frame[Y0:Y1, X0:X1].astype(np.float32)

				with self.lock:
					self.frame = roi

	def get(self):
		with self.lock:
			return None if self.frame is None else self.frame.copy()

class Server:
	def __init__(self, port=54321):
		self.port = port
		self.cam = Camera()
		self.idx = 0

	def start(self):
		srv = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
		srv.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
		srv.bind(("0.0.0.0", self.port))
		srv.listen(10)

		print("Waiting for clients...")

		while True:
			conn, addr = srv.accept()
			print("[CONNECTION]", addr)

			threading.Thread(
				target=self._client,
				args=(conn, addr),
				daemon=True
			).start()

	def _client(self, conn, addr):
		try:
			while True:
				frame = self.cam.get()
				if frame is None:
					frame = TV_BARS

				header = struct.pack(
					_CAM_HDR_FMT,
					_CAM_MAGIC,
					H,
					W,
					self.idx,
					time.time()
				)

				conn.sendall(header)
				conn.sendall(frame.tobytes())

				self.idx += 1
				time.sleep(1 / 25)

		except:
			pass
		finally:
			conn.close()
			print("[DISCONNECTION]", addr)

if __name__ == "__main__":
	Server().start()
