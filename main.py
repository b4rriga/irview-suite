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

from __future__ import annotations
import ctypes, json, os, socket, struct, sys, threading, time, warnings
from typing import Optional

import numpy as np
from scipy.linalg import svd as scipy_svd

import matplotlib
matplotlib.use("Qt5Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.colors import LinearSegmentedColormap

from PyQt5.QtWidgets import (
	QApplication, QButtonGroup, QCheckBox, QComboBox, QDialog,
	QDialogButtonBox, QFileDialog, QFrame, QGridLayout, QGroupBox,
	QHBoxLayout, QLabel, QLineEdit, QMainWindow, QMessageBox,
	QPushButton, QScrollArea, QShortcut, QSizePolicy, QSlider,
	QSpinBox, QDoubleSpinBox, QSplitter, QVBoxLayout, QWidget,
)
from PyQt5.QtCore import Qt, QTimer, pyqtSignal
from PyQt5.QtGui import (QDoubleValidator, QImage, QKeySequence,
                         QPainter, QPixmap, QColor, QPen, QFont)

# ---------------------------------------------------------------------------
# C extension
# ---------------------------------------------------------------------------

def _load_ircore():
	so = os.path.join(os.path.dirname(os.path.abspath(__file__)), "ircore.so")
	if not os.path.exists(so): return None
	try:
		lib = ctypes.CDLL(so)
		fp  = ctypes.POINTER(ctypes.c_float)
		lib.hmf.argtypes             = [fp, fp, ctypes.c_int, ctypes.c_int]
		lib.gauss.argtypes           = [fp, fp, ctypes.c_int, ctypes.c_int, ctypes.c_float]
		lib.corrct.argtypes          = [fp, fp, ctypes.c_int, ctypes.c_int, ctypes.c_int]
		lib.pct_standardize.argtypes = [fp, ctypes.c_int, ctypes.c_int]
		lib.extrap_sib.argtypes      = [fp, ctypes.c_float, ctypes.c_float, fp, ctypes.c_int]
		lib.minmax.argtypes          = [fp, ctypes.c_int,
											ctypes.POINTER(ctypes.c_float),
											ctypes.POINTER(ctypes.c_float)]
		for fn in (lib.hmf, lib.gauss, lib.corrct,
				   lib.pct_standardize, lib.extrap_sib, lib.minmax):
			fn.restype = None
		return lib
	except Exception: return None

_lib = _load_ircore()
_fp  = ctypes.POINTER(ctypes.c_float)
def _ptr(a): return a.ctypes.data_as(_fp)

def _nanminmax(arr):
	f32 = np.ascontiguousarray(arr.ravel(), dtype=np.float32)
	if _lib:
		lo = ctypes.c_float(); hi = ctypes.c_float()
		_lib.minmax(_ptr(f32), len(f32), ctypes.byref(lo), ctypes.byref(hi))
		return float(lo.value), float(hi.value)
	return float(np.nanmin(arr)), float(np.nanmax(arr))

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_VJ = np.array([1/12,-385/12,1279,-46871/3,505465/6,
				-473915/2,1127735/3,-1020215/3,328125/2,-65625/2])
_MAGIC       = b"IRVW"
_CAM_MAGIC   = b"IRFR"
_CAM_HDR_FMT = "<4sIIId"
_CAM_HDR_LEN = struct.calcsize(_CAM_HDR_FMT)

CLAMP_IMG, CLAMP_SEQ, CLAMP_ROI, CLAMP_MAN = 0, 1, 2, 3
MODE_IDLE, MODE_FILE, MODE_STREAM = 0, 1, 2

# ---------------------------------------------------------------------------
# Colourmaps
# ---------------------------------------------------------------------------

def _make_iron():
	stops = [(0.,(0,0,0)),(0.25,(.5,0,.5)),(.5,(1,0,0)),(.75,(1,.6,0)),(1.,(1,1,1))]
	xs = [s[0] for s in stops]; cs = [s[1] for s in stops]
	return LinearSegmentedColormap('iron',
		{ch:[(x,c[i],c[i]) for x,c in zip(xs,cs)]
		 for i,ch in enumerate(('red','green','blue'))})

plt.colormaps.register(_make_iron(), force=True)

_CMAPS = [
	("autumn","Autumn",256), ("bone","Bone",256), ("gist_ncar","Colorcube",256),
	("cool","Cool",256),     ("copper","Copper",256), ("gray","Gray",256),
	("hot","Hot",256),       ("hsv","HSV",256),   ("iron","Iron",256),
	("jet","Jet",256),       ("pink","Pink",256), ("spring","Spring",256),
	("summer","Summer",256), ("winter","Winter",256),
]
_CMAP_KEYS    = [c[0] for c in _CMAPS]
_CMAP_LABELS  = [c[1] for c in _CMAPS]
_CMAP_MAXCOLS = [c[2] for c in _CMAPS]
_EXTRAP_TYPES = ["Semi Infinite Body","Slab (adiabatic)","Slab (with losses)"]

# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------

def _ensure3d(a): return a[:,:,np.newaxis] if a.ndim==2 else a

def load_mat(path):
	import scipy.io
	mat = scipy.io.loadmat(path)
	cands = {k:v for k,v in mat.items()
			 if not k.startswith('_') and isinstance(v,np.ndarray) and v.ndim>=2}
	if not cands: raise ValueError("No 2-D+ array in .mat file.")
	arr = np.asarray(max(cands.values(), key=lambda a:a.size), dtype=np.float32)
	if arr.ndim==3 and arr.shape[0]<arr.shape[1] and arr.shape[0]<arr.shape[2]:
		arr = np.moveaxis(arr, 0, -1)
	return _ensure3d(arr)

def save_frame_mat(path, frame):
	import scipy.io
	scipy.io.savemat(path, {"frame": np.asarray(frame, dtype=np.float64)})

def save_capture_mat(path, cube):
	"""cube: (H, W, N) float32"""
	import scipy.io
	scipy.io.savemat(path, {"capture": np.asarray(cube, dtype=np.float64)})

def save_irv(path, state, seq):
	payload = json.dumps(state).encode()
	arr = np.ascontiguousarray(seq, dtype=np.float32); H,W,N = arr.shape
	with open(path,"wb") as f:
		f.write(_MAGIC); f.write(struct.pack("<I",len(payload))); f.write(payload)
		f.write(struct.pack("<III",H,W,N)); f.write(arr.tobytes())

def load_irv(path):
	with open(path,"rb") as f:
		if f.read(4)!=_MAGIC: raise ValueError("Not a valid .irv file.")
		slen = struct.unpack("<I",f.read(4))[0]
		state = json.loads(f.read(slen)); rest = f.read()
	H,W,N = struct.unpack("<III",rest[:12])
	return state, np.frombuffer(rest[12:],dtype=np.float32).reshape(H,W,N).copy()

# ---------------------------------------------------------------------------
# Algorithms
# ---------------------------------------------------------------------------

def smthgauss(img, v):
	f32 = np.ascontiguousarray(img, dtype=np.float32)
	if _lib:
		dst = np.empty_like(f32)
		_lib.gauss(_ptr(f32),_ptr(dst),f32.shape[0],f32.shape[1],ctypes.c_float(v))
		return dst.astype(img.dtype)
	from scipy.signal import convolve2d
	i_max=10*np.sqrt(v); i=np.arange(0,round(i_max)+1,dtype=float)
	y=np.exp(-((i-i_max/2)**2)/(2*v))/np.sqrt(2*np.pi*v); y/=y.sum()
	bw=len(y)//2; pad=np.pad(img,bw,mode="edge")
	return convolve2d(pad,np.outer(y,y),mode="same")[bw:-bw,bw:-bw]

def hmf(A):
	f32 = np.ascontiguousarray(A, dtype=np.float32)
	if _lib:
		dst = np.empty_like(f32)
		_lib.hmf(_ptr(f32),_ptr(dst),f32.shape[0],f32.shape[1])
		return dst.astype(A.dtype)
	def _g(arr,offs):
		h,w=arr.shape; rs,cs=np.mgrid[0:h,0:w]; rs,cs=rs.ravel(),cs.ravel()
		out=np.empty((len(rs),len(offs)),dtype=arr.dtype)
		for k,(dr,dc) in enumerate(offs):
			out[:,k]=arr[np.clip(rs+dr,0,h-1),np.clip(cs+dc,0,w-1)]
		return out
	cr=[(-2,0),(-1,0),(0,-2),(0,-1),(0,0),(0,1),(0,2),(1,0),(2,0)]
	di=[(-2,-2),(-1,-1),(0,0),(1,1),(2,2),(-2,2),(-1,1),(1,-1),(2,-2)]
	h,w=A.shape; Af=A.astype(float)
	return np.median(np.stack([np.median(_g(Af,cr),axis=1),
								np.median(_g(Af,di),axis=1),
								Af.ravel()],axis=1),axis=1).reshape(h,w).astype(A.dtype)

def correlated_contrast(seq):
	try:
		H,W,N = seq.shape
		f32 = np.ascontiguousarray(seq.reshape(H*W,N), dtype=np.float32)
		if _lib:
			out = np.empty(H*W,dtype=np.float32)
			_lib.corrct(_ptr(f32),_ptr(out),H,W,N)
			return out.reshape(H,W)
		sf=seq.astype(float); ref=sf.mean(axis=(0,1)); CC=np.zeros((H,W))
		for i in range(H):
			p=sf[i,:,:].T; pc=p-p.mean(0); rc=ref-ref.mean()
			den=np.sqrt((pc**2).sum(0))*np.sqrt((rc**2).sum()); den[den==0]=1.
			CC[i,:]=1.-(pc*rc[:,None]).sum(0)/den
		return np.cbrt(CC**5).astype(np.float32)
	except MemoryError: return None

def pct_standardize(m):
	f32 = np.ascontiguousarray(m, dtype=np.float32)
	if _lib:
		_lib.pct_standardize(_ptr(f32), f32.shape[0], f32.shape[1])
		return f32
	mu=m.mean(0); s=m.std(0,ddof=1); s[s==0]=1.
	return ((m-mu)/s).astype(np.float32)

def _sp(t):
	t=np.atleast_1d(t).ravel()
	return np.arange(1,len(_VJ)+1,dtype=float).reshape(-1,1)*np.log(2)/t.reshape(1,-1)

def _sFt(Fp,t):
	return np.log(2)/np.atleast_1d(t).ravel()*np.sum(_VJ.reshape(-1,1)*Fp,axis=0)

def extrapolated_temperature(t_p,T_p,xdata,etype,L=4e-3,diff=0.1e-6,cond=0.2,h=10.):
	t=np.atleast_1d(np.asarray(xdata,float)).ravel()
	T=np.asarray(T_p,float); sc=t.size==1
	with warnings.catch_warnings():
		warnings.simplefilter("ignore",RuntimeWarning)
		if etype==1:
			if sc and _lib and T.ndim==2:
				f32=np.ascontiguousarray(T,dtype=np.float32); dst=np.empty_like(f32)
				_lib.extrap_sib(_ptr(f32),ctypes.c_float(t_p),ctypes.c_float(float(t[0])),
									_ptr(dst),f32.size)
				return dst
			r=np.sqrt(t_p)/np.sqrt(t); return T*float(r[0]) if sc else T*r
		def _q2(tv): p=_sp(tv); return _sFt((1/np.sqrt(p))/np.tanh(np.sqrt(p*L**2/diff)),tv)
		def _q3(tv):
			p=_sp(tv); k=np.sqrt(p/diff)
			return _sFt((h*np.sinh(k*L)/(cond*k)+np.cosh(k*L))/(
				2*np.cosh(k*L)*h+cond*k*np.sinh(k*L)+h**2*np.sinh(k*L)/(cond*k)),tv)
		fn=_q2 if etype==2 else _q3
		r=fn(t)/fn(np.array([t_p]))[0]; return T*float(r[0]) if sc else T*r

def apply_filter(img, ctrl):
	if ctrl.btn_median.isChecked(): return hmf(img)
	if ctrl.btn_gauss.isChecked() and ctrl.variance>0: return smthgauss(img,ctrl.variance)
	return img

# ---------------------------------------------------------------------------
# Higher-level processing algorithms (MATLAB script ports)
# ---------------------------------------------------------------------------

def absolute_contrast(seq, cold, ref_row, ref_col):
	"""
	CA: Absolute Contrast.
	seq  : (H,W,N) float32
	cold : (H,W)   float32  cold (pre-flash) frame
	Returns (H,W,N) float32.
	"""
	H,W,N = seq.shape
	f32s = np.ascontiguousarray(seq.reshape(H*W,N), dtype=np.float32)
	f32c = np.ascontiguousarray(cold.ravel(),        dtype=np.float32)
	ref_px = int(ref_row)*W + int(ref_col)
	out  = np.empty_like(f32s)
	if _lib:
		_lib.ca.argtypes = [_fp,_fp,ctypes.c_int,_fp,ctypes.c_int,ctypes.c_int]
		_lib.ca.restype  = None
		_lib.ca(_ptr(f32s),_ptr(f32c),ctypes.c_int(ref_px),_ptr(out),
					ctypes.c_int(H*W),ctypes.c_int(N))
	else:
		B = seq.astype(float) - cold[:,:,np.newaxis]
		ref_sig = B[ref_row,ref_col,:]
		out_np  = B - ref_sig[np.newaxis,np.newaxis,:]
		return out_np.astype(np.float32)
	return out.reshape(H,W,N)

def differential_absolute_contrast(seq, cold, t_prime):
	"""
	DAC: Differential Absolute Contrast.
	t_prime : 1-based reference frame index (MATLAB convention).
	Returns (H,W, N-t_prime) float32.
	"""
	H,W,N = seq.shape
	out_N  = N - t_prime
	if out_N <= 0: raise ValueError(f"t_prime={t_prime} >= N={N}")
	f32s = np.ascontiguousarray(seq.reshape(H*W,N), dtype=np.float32)
	f32c = np.ascontiguousarray(cold.ravel(),        dtype=np.float32)
	out  = np.empty((H*W,out_N), dtype=np.float32)
	if _lib:
		_lib.dac.argtypes = [_fp,_fp,ctypes.c_int,_fp,ctypes.c_int,ctypes.c_int]
		_lib.dac.restype  = None
		_lib.dac(_ptr(f32s),_ptr(f32c),ctypes.c_int(t_prime),
					 _ptr(out),ctypes.c_int(H*W),ctypes.c_int(N))
	else:
		B = seq.astype(float) - cold[:,:,np.newaxis]
		tp0 = t_prime-1
		for i in range(out_N):
			out_np = B[:,:,i+tp0] - np.sqrt(t_prime/(i+1))*B[:,:,tp0]
			out[: ,i] = out_np.ravel()
	return out.reshape(H,W,out_N)

def rx_detector(seq):
	"""
	RX anomaly detector.  Returns (H,W) float32 score map.
	Spheres the data (eigen-whitening) then computes ||Mx||^2 per pixel.
	"""
	H,W,N = seq.shape
	m  = seq.reshape(H*W,N).astype(float)
	m -= m.mean(axis=1, keepdims=True)
	# Covariance of temporal signals
	cov = np.cov(m, rowvar=False)               # (N,N)
	vals,vecs = np.linalg.eigh(cov)
	vals = np.maximum(vals, 1e-12)
	sphere = (vecs / np.sqrt(vals)).T           # (N,N) sphering matrix
	m  = np.ascontiguousarray(m,      dtype=np.float32)
	sp = np.ascontiguousarray(sphere,  dtype=np.float32)
	out    = np.empty(H*W, dtype=np.float32)
	if _lib:
		_lib.rx.argtypes=[_fp,_fp,_fp,ctypes.c_int,ctypes.c_int]
		_lib.rx.restype=None
		_lib.rx(_ptr(m),_ptr(sp),_ptr(out),ctypes.c_int(H*W),ctypes.c_int(N))
	else:
		sphered = (sphere @ m.T).T              # (H*W, N)
		out_np  = (sphered**2).sum(axis=1)
		out     = out_np.astype(np.float32)
	return out.reshape(H,W)

def skq_stats(seq):
	"""
	SKQ: per-pixel kurtosis, skewness, RX (un-sphered), 5th moment.
	Returns dict with keys 'kurtosis','skewness','rx','quinto', each (H,W) float32.
	"""
	H,W,N = seq.shape
	m  = seq.reshape(H*W,N).astype(float)
	m -= m.mean(axis=1, keepdims=True)
	f32 = np.ascontiguousarray(m, dtype=np.float32)
	kurt  = np.empty(H*W,dtype=np.float32)
	skew  = np.empty(H*W,dtype=np.float32)
	rxout = np.empty(H*W,dtype=np.float32)
	quint = np.empty(H*W,dtype=np.float32)
	if _lib:
		_lib.skq.argtypes=[_fp,_fp,_fp,_fp,_fp,ctypes.c_int,ctypes.c_int]
		_lib.skq.restype=None
		_lib.skq(_ptr(f32),_ptr(kurt),_ptr(skew),_ptr(rxout),_ptr(quint),
					 ctypes.c_int(H*W),ctypes.c_int(N))
	else:
		sig2 = m.var(axis=1);  sig2[sig2<1e-30]=1e-30
		sig  = np.sqrt(sig2)
		kurt[:]  = (m**4).mean(axis=1)/sig2**2
		skew[:]  = (m**3).mean(axis=1)/sig2**1.5
		rxout[:] = (m**2).sum(axis=1)
		quint[:] = (m**5).mean(axis=1)/(sig2**2*sig)
	return {"kurtosis":kurt.reshape(H,W),"skewness":skew.reshape(H,W),
			"rx":rxout.reshape(H,W),"quinto":quint.reshape(H,W)}

def pca_mmethod(seq):
	"""
	PCA M-Method: SVD, keep components until cumulative variance > 99.9%.
	Returns (U, s_norm, n_components) where U is (H,W, n_components).
	"""
	H,W,N = seq.shape
	m = pct_standardize(np.ascontiguousarray(seq.reshape(H*W,N), dtype=np.float32))
	U,sv,_ = scipy_svd(m.astype(float), full_matrices=False)
	cmax = sv.sum(); cumsum = 0.0; I = len(sv)
	for i,v in enumerate(sv):
		cumsum += v
		if cumsum/cmax > 0.999: I=i+1; break
	return U[:,:I].reshape(H,W,I), sv/sv.sum(), I

def tsr_polyfit(seq, cold, degree=5):
	"""
	TSR: fit log(T-T_cold) vs log(t) with a polynomial of given degree.
	Returns coefficient cube (H,W, degree+1) float32 and synthetic data (H,W,N) float32.
	"""
	H,W,N = seq.shape
	t_log = np.log(np.arange(1,N+1,dtype=float))
	B     = np.log(np.maximum(seq.astype(float)-cold[:,:,np.newaxis], 1e-10))
	# Vandermonde matrix (N, degree+1)
	V  = np.vander(t_log, degree+1, increasing=True)   # (N, d+1)
	# Solve per pixel via lstsq on the reshaped matrix
	B2 = B.reshape(H*W, N)                             # (H*W, N)
	# coeffs: (H*W, degree+1)
	coeffs, _,_,_ = np.linalg.lstsq(V, B2.T, rcond=None)  # (d+1, H*W)
	coeffs = coeffs.T.astype(np.float32)                    # (H*W, d+1)
	synth  = (V @ coeffs.T).T.reshape(H,W,N).astype(np.float32)
	return coeffs.reshape(H,W,degree+1), synth

def haar_dwt(seq):
	"""
	1-level Haar DWT per pixel along time axis.
	Returns (H,W, N//2) approximation and (H,W, N//2) detail coefficients.
	"""
	H,W,N = seq.shape
	if N%2!=0: seq=seq[:,:,:N-1]; N=N-1
	f32 = np.ascontiguousarray(seq.reshape(H*W,N), dtype=np.float32)
	ca  = np.empty((H*W,N//2),dtype=np.float32)
	cd  = np.empty((H*W,N//2),dtype=np.float32)
	if _lib:
		_lib.haar_dwt.argtypes=[_fp,_fp,_fp,ctypes.c_int,ctypes.c_int]
		_lib.haar_dwt.restype=None
		_lib.haar_dwt(_ptr(f32),_ptr(ca),_ptr(cd),ctypes.c_int(H*W),ctypes.c_int(N))
	else:
		S  = 0.7071067811865476
		ca = S*(seq[:,:,0::2]+seq[:,:,1::2]).astype(np.float32).reshape(H*W,N//2)
		cd = S*(seq[:,:,0::2]-seq[:,:,1::2]).astype(np.float32).reshape(H*W,N//2)
	return ca.reshape(H,W,N//2), cd.reshape(H,W,N//2)

# ---------------------------------------------------------------------------
# Colourmap lookup table (numpy, used by OpenCV canvas)
# ---------------------------------------------------------------------------

def _build_lut(cmap_name, n_colors):
	"""Return uint8 (256,3) BGR lookup table for the given matplotlib cmap."""
	cmap = plt.get_cmap(cmap_name).resampled(n_colors)
	idx  = np.linspace(0, 1, 256)
	rgba = (cmap(idx)[:,:3]*255).astype(np.uint8)  # (256,3) RGB
	return rgba[:, ::-1]                            # -> BGR for OpenCV / QImage path

# ---------------------------------------------------------------------------
# Camera receiver (TCP)
# ---------------------------------------------------------------------------

class CameraReceiver:
	def __init__(self):
		self._lock    = threading.Lock()
		self._latest  : Optional[np.ndarray] = None
		self._ts      = 0.0; self._idx = 0
		self._sock    : Optional[socket.socket] = None
		self._thread  : Optional[threading.Thread] = None
		self._running = False
		self.H=self.W=0; self.error=""

	def start(self, host, port):
		self.stop(); self._running=True; self.error=""
		self._thread = threading.Thread(target=self._loop, args=(host,port), daemon=True)
		self._thread.start()

	def stop(self):
		self._running=False
		if self._sock:
			try: self._sock.close()
			except: pass
			self._sock=None
		if self._thread: self._thread.join(timeout=2.0); self._thread=None
		with self._lock: self._latest=None

	def is_running(self): return self._running

	def get_latest(self):
		with self._lock:
			return (self._latest.copy() if self._latest is not None else None,
					self._idx, self._ts)

	def _loop(self, host, port):
		try:
			self._sock=socket.socket(socket.AF_INET, socket.SOCK_STREAM)
			self._sock.connect((host,int(port))); self._sock.settimeout(5.0)
		except Exception as e:
			self.error=str(e); self._running=False; return
		while self._running:
			try:
				hdr=self._recvall(_CAM_HDR_LEN)
				if hdr is None: break
				magic,H,W,fidx,ts=struct.unpack(_CAM_HDR_FMT,hdr)
				if magic!=_CAM_MAGIC: break
				payload=self._recvall(H*W*4)
				if payload is None: break
				frame=np.frombuffer(payload,dtype=np.float32).reshape(H,W).copy()
				with self._lock:
					self._latest=frame; self._idx=fidx; self._ts=ts
					self.H=H; self.W=W
			except socket.timeout: continue
			except Exception: break
		self._running=False

	def _recvall(self,n):
		buf=bytearray()
		while len(buf)<n:
			try: chunk=self._sock.recv(n-len(buf))
			except: return None
			if not chunk: return None
			buf.extend(chunk)
		return bytes(buf)

# ---------------------------------------------------------------------------
# Stylesheet
# ---------------------------------------------------------------------------

_STYLE = """
QWidget {
	background-color: #e4e4e4;
	color: #1a1a1a;
	font-size: 13px;
}
QGroupBox {
	background-color: #dedede;
	border: 1px solid #a8a8a8;
	border-radius: 4px;
	margin-top: 16px;
	padding: 10px 6px 6px 6px;
	font-weight: bold;
	font-size: 13px;
}
QGroupBox::title {
	subcontrol-origin: margin;
	subcontrol-position: top left;
	left: 10px;
	padding: 0 5px;
	color: #2a2a2a;
	background: transparent;
}
QPushButton {
	background-color: #d4d4d4;
	border: 1px solid #979797;
	border-radius: 3px;
	padding: 5px 11px;
	min-height: 28px;
	color: #1a1a1a;
	font-size: 13px;
}
QPushButton:hover    { background-color: #c6c6c6; border-color: #707070; }
QPushButton:pressed  { background-color: #b8b8b8; }
QPushButton:checked  { background-color: #4a86c8; border-color: #2a66a8; color: #ffffff; }
QPushButton:disabled { background-color: #d2d2d2; border-color: #b8b8b8; color: #909090; }
QSlider::groove:horizontal {
	height: 5px; background: #b8b8b8; border-radius: 2px; margin: 0 3px;
}
QSlider::groove:horizontal:disabled { background: #cccccc; }
QSlider::handle:horizontal {
	width: 15px; height: 15px; margin: -5px 0;
	background: #6699cc; border: 1px solid #4477aa; border-radius: 8px;
}
QSlider::handle:horizontal:disabled { background: #bbbbbb; border-color: #aaaaaa; }
QCheckBox { spacing: 8px; min-height: 24px; color: #1a1a1a; font-size: 13px; }
QCheckBox::indicator {
	width: 15px; height: 15px;
	border: 1px solid #979797; border-radius: 2px; background: #f0f0f0;
}
QCheckBox::indicator:checked  { background: #4a86c8; border-color: #2a66a8; }
QCheckBox::indicator:disabled { background: #d8d8d8; border-color: #bbbbbb; }
QLabel { color: #1a1a1a; padding: 2px 0; background: transparent; font-size: 13px; }
QComboBox {
	background: #efefef; border: 1px solid #979797; border-radius: 3px;
	padding: 4px 24px 4px 8px; min-height: 28px; color: #1a1a1a; font-size: 13px;
}
QComboBox::drop-down {
	subcontrol-origin: padding; subcontrol-position: top right;
	width: 22px; border-left: 1px solid #b0b0b0;
}
QComboBox::down-arrow {
	image: none;
	border-left: 5px solid transparent;
	border-right: 5px solid transparent;
	border-top: 6px solid #555555;
	width: 0; height: 0;
}
QComboBox:disabled { background: #e0e0e0; color: #888888; }
QLineEdit {
	background: #f8f8f8; border: 1px solid #979797; border-radius: 3px;
	padding: 4px 7px; min-height: 28px; color: #1a1a1a; font-size: 13px;
}
QLineEdit:disabled { background: #e8e8e8; color: #888888; }
QSpinBox, QDoubleSpinBox {
	background: #f8f8f8; border: 1px solid #979797; border-radius: 3px;
	padding: 4px 7px; min-height: 28px; color: #1a1a1a; font-size: 13px;
}
QScrollArea { border: none; background: transparent; }
QScrollBar:vertical { background: #d8d8d8; width: 12px; border: none; }
QScrollBar::handle:vertical {
	background: #aaaaaa; border-radius: 5px; min-height: 20px;
}
QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical { height: 0; }
QMessageBox { font-size: 13px; }
QInputDialog { font-size: 13px; }
QInputDialog QLabel { min-width: 260px; padding: 8px 4px; }
QInputDialog QLineEdit { min-width: 260px; }
QFrame[frameShape="4"], QFrame[frameShape="5"] { color: #b0b0b0; }
"""

# ---------------------------------------------------------------------------
# Fast OpenCV-style main display canvas (QLabel + QPixmap, no Matplotlib)
# ---------------------------------------------------------------------------

class IRDisplay(QWidget):
	"""
	High-performance IR image display using QPixmap rendering.
	Supports: colour mapping, axis ticks with real coordinates,
	grid overlay, x/y inversion, ROI rubber-band overlay.
	Mouse signals mirror the old IRCanvas API.
	"""
	canvas_clicked   = pyqtSignal(float, float)   # display-space coords
	canvas_mouse_pos = pyqtSignal(float, float)
	canvas_mouse_out = pyqtSignal()

	def __init__(self):
		super().__init__()
		self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
		self.setMouseTracking(True)
		self.setMinimumSize(200, 150)

		self._img    : Optional[np.ndarray] = None   # (H,W) float32 current frame
		self._lut    : Optional[np.ndarray] = None   # (256,3) uint8 BGR
		self._lo     = 0.0
		self._hi     = 1.0
		self._inv_x  = False
		self._inv_y  = False
		self._show_axis = True
		self._show_grid = False
		# ROI offset: (x0, y0) in full-frame pixels, or None
		self._roi_offset: Optional[tuple] = None
		# ROI rubber-band overlay (display-space)
		self._rb_active = False
		self._rb_pt1: Optional[tuple] = None
		self._rb_pt2: Optional[tuple] = None
		self._rb_shading = True

		# Margins for axis ticks (pixels)
		self._MARG_L = 54
		self._MARG_T = 30
		self._MARG_R = 14
		self._MARG_B = 24

		# Colourbar
		self._CBAR_W = 26
		self._CBAR_LABEL_W = 58

	# -- public API -----------------------------------------------------------

	def show_frame(self, img, clim, lut, roi_offset=None):
		"""img: (H,W) float, clim: (lo,hi), lut: (256,3) uint8 BGR."""
		self._img = img
		self._lo, self._hi = clim
		self._lut = lut
		self._roi_offset = roi_offset
		self.update()

	def set_clim(self, lo, hi):
		self._lo = lo; self._hi = hi; self.update()

	def set_invert_x(self, on): self._inv_x = on; self.update()
	def set_invert_y(self, on): self._inv_y = on; self.update()
	def toggle_axis(self, on): self._show_axis = on; self.update()
	def toggle_grid(self, on): self._show_grid = on; self.update()

	def roi_show_pt1(self, x, y):
		self._rb_active = True; self._rb_pt1 = (x,y); self._rb_pt2 = None
		self.update()

	def roi_update_preview(self, x0, y0, x1, y1, img_h, img_w,
						   x_offset=0, y_offset=0):
		self._rb_active = True
		self._rb_pt1 = (x0, y0); self._rb_pt2 = (x1, y1)
		self.update()

	def roi_clear_overlay(self):
		self._rb_active = False; self._rb_pt1 = self._rb_pt2 = None
		self.update()

	# -- coordinate transforms ------------------------------------------------

	def _img_rect(self):
		"""Return (x, y, w, h) of the image area inside the widget."""
		ml = self._MARG_L if self._show_axis else 6
		mt = self._MARG_T if self._show_axis else 6
		mr = self._MARG_R
		mb = self._MARG_B if self._show_axis else 6
		cb = self._CBAR_W + self._CBAR_LABEL_W + 6
		W = self.width()  - ml - mr - cb
		H = self.height() - mt - mb
		return ml, mt, max(W, 1), max(H, 1)

	def _display_to_img(self, dx, dy):
		"""Convert widget pixel (dx,dy) -> image float coords."""
		ml, mt, iw, ih = self._img_rect()
		if self._img is None: return dx, dy
		H, W = self._img.shape
		fx = (dx - ml) / iw * W
		fy = (dy - mt) / ih * H
		if self._inv_x: fx = W - 1 - fx
		if self._inv_y: fy = H - 1 - fy
		return fx, fy

	# -- painting --------------------------------------------------------------

	def paintEvent(self, ev):
		p = QPainter(self)
		p.setRenderHint(QPainter.Antialiasing, False)
		ml, mt, iw, ih = self._img_rect()

		# Background
		p.fillRect(self.rect(), QColor("#e4e4e4"))

		if self._img is None or self._lut is None:
			p.end(); return

		H, W = self._img.shape

		# -- remap float32 -> uint8 via LUT -----------------------------------
		span = self._hi - self._lo
		if span == 0: span = 1e-9
		normed = np.clip((self._img - self._lo) / span * 255, 0, 255).astype(np.uint8)
		if self._inv_x: normed = normed[:, ::-1]
		if self._inv_y: normed = normed[::-1, :]
		rgb = self._lut[normed]             # (H, W, 3) uint8 BGR
		rgb_c = np.ascontiguousarray(rgb[:,:,::-1])   # -> RGB
		qi = QImage(rgb_c.tobytes(), W, H, W*3, QImage.Format_RGB888)
		px = QPixmap.fromImage(qi).scaled(iw, ih, Qt.IgnoreAspectRatio,
										  Qt.FastTransformation)
		p.drawPixmap(ml, mt, px)

		# -- ROI rubber-band overlay ----------------------------------------
		if self._rb_active and self._rb_pt1 is not None:
			x0, y0 = self._rb_pt1
			if self._rb_pt2 is None:
				# Just show first point marker
				px0 = int(ml + x0/W*iw); py0 = int(mt + y0/H*ih)
				p.setPen(QPen(QColor(0,255,0), 2))
				p.setBrush(QColor(0,255,0,180))
				p.drawEllipse(px0-5, py0-5, 10, 10)
			else:
				x1, y1 = self._rb_pt2
				# Clamp to image space
				rx0 = int(ml + np.clip(min(x0,x1),0,W-1)/W*iw)
				rx1 = int(ml + np.clip(max(x0,x1),0,W-1)/W*iw)
				ry0 = int(mt + np.clip(min(y0,y1),0,H-1)/H*ih)
				ry1 = int(mt + np.clip(max(y0,y1),0,H-1)/H*ih)
				# Shading outside selection
				shade = QColor(0,0,0,100)
				for (sx,sy,sw,sh) in [
					(ml,       mt,           iw,        ry0-mt),
					(ml,       ry1,          iw,        mt+ih-ry1),
					(ml,       ry0,          rx0-ml,    ry1-ry0),
					(rx1,      ry0,          ml+iw-rx1, ry1-ry0),
				]:
					if sw>0 and sh>0:
						p.fillRect(sx,sy,sw,sh,shade)
				# Lime rectangle
				p.setPen(QPen(QColor(0,255,0), 2))
				p.setBrush(Qt.NoBrush)
				p.drawRect(rx0, ry0, rx1-rx0, ry1-ry0)

		# -- axis ticks and grid ----------------------------------------------
		N_TICKS = 5
		ox = (self._roi_offset[0] if self._roi_offset else 0)
		oy = (self._roi_offset[1] if self._roi_offset else 0)
		# Precompute tick fractions and pixel positions
		x_fracs = [i / (N_TICKS - 1) for i in range(N_TICKS)]
		y_fracs = [i / (N_TICKS - 1) for i in range(N_TICKS)]
		x_wx    = [int(ml + f * iw) for f in x_fracs]
		y_wy    = [int(mt + f * ih) for f in y_fracs]

		if self._show_grid:
			pen = QPen(QColor(80, 80, 80, 80)); pen.setWidth(1)
			p.setPen(pen)
			for wx in x_wx[1:-1]:
				p.drawLine(wx, mt, wx, mt + ih)
			for wy in y_wy[1:-1]:
				p.drawLine(ml, wy, ml + iw, wy)

		if self._show_axis:
			p.setPen(QPen(QColor("#1a1a1a"), 1))
			font = QFont(); font.setPixelSize(12); p.setFont(font)

			# X axis: top when inv_y (row 0 at top), bottom otherwise
			x_on_top = self._inv_y
			x_tick_y = mt if x_on_top else mt + ih
			x_label_y = (x_tick_y - 22) if x_on_top else (x_tick_y + 4)
			for frac, wx in zip(x_fracs, x_wx):
				img_col = frac * (W - 1) if not self._inv_x else (1 - frac) * (W - 1)
				real_x  = int(round(img_col)) + ox
				if x_on_top:
					p.drawLine(wx, x_tick_y - 6, wx, x_tick_y)
				else:
					p.drawLine(wx, x_tick_y, wx, x_tick_y + 6)
				p.drawText(wx - 22, x_label_y, 44, 18, Qt.AlignCenter, str(real_x))

			# Y axis: left when not inv_x, right when inv_x
			y_on_left = not self._inv_x
			y_tick_x  = ml if y_on_left else ml + iw
			if y_on_left:
				label_rect = lambda wy: (0, wy - 10, ml - 6, 20)
				label_align = Qt.AlignRight | Qt.AlignVCenter
			else:
				label_rect = lambda wy: (ml + iw + 6, wy - 10, self._MARG_R + 6, 20)
				label_align = Qt.AlignLeft | Qt.AlignVCenter
			for frac, wy in zip(y_fracs, y_wy):
				img_row = frac * (H - 1) if not self._inv_y else (1 - frac) * (H - 1)
				real_y  = int(round(img_row)) + oy
				if y_on_left:
					p.drawLine(y_tick_x - 6, wy, y_tick_x, wy)
				else:
					p.drawLine(y_tick_x, wy, y_tick_x + 6, wy)
				rx, ry, rw, rh = label_rect(wy)
				p.drawText(rx, ry, rw, rh, label_align, str(real_y))

			# Border
			p.setPen(QPen(QColor("#777777"), 1))
			p.setBrush(Qt.NoBrush)
			p.drawRect(ml, mt, iw, ih)

		# -- colorbar ---------------------------------------------------------
		cb_x = self.width() - self._CBAR_W - self._CBAR_LABEL_W - 4
		cb_h = ih
		cb_idx = np.arange(255, -1, -1, dtype=np.uint8)          # (256,) 1D indices
		cb_rgb = np.ascontiguousarray(self._lut[cb_idx][:,::-1])  # (256,3) RGB
		# Each row is one colour; repeat to width=1 px, then scale
		cb_row = cb_rgb.reshape(256, 1, 3)
		cb_c   = np.ascontiguousarray(np.repeat(cb_row, 1, axis=1))
		qi_cb  = QImage(cb_c.tobytes(), 1, 256, 3, QImage.Format_RGB888)
		px_cb  = QPixmap.fromImage(qi_cb).scaled(
			self._CBAR_W, cb_h, Qt.IgnoreAspectRatio, Qt.FastTransformation)
		p.drawPixmap(cb_x, mt, px_cb)
		p.setPen(QPen(QColor("#777777"),1))
		p.drawRect(cb_x, mt, self._CBAR_W, cb_h)

		# Colorbar labels (7 evenly spaced)
		font2 = QFont(); font2.setPixelSize(11); p.setFont(font2)
		p.setPen(QPen(QColor("#1a1a1a"),1))
		for i in range(7):
			frac = i / 6
			val  = self._hi + (self._lo - self._hi) * frac
			wy   = int(mt + frac * cb_h)
			p.drawText(cb_x + self._CBAR_W + 4, wy - 9, self._CBAR_LABEL_W - 4, 18,
					   Qt.AlignLeft | Qt.AlignVCenter, f"{val:.3g}")

		p.end()

	# -- mouse events ----------------------------------------------------------

	def mousePressEvent(self, ev):
		if ev.button() == Qt.LeftButton:
			fx, fy = self._display_to_img(ev.x(), ev.y())
			self.canvas_clicked.emit(fx, fy)

	def mouseMoveEvent(self, ev):
		ml, mt, iw, ih = self._img_rect()
		if (ml <= ev.x() <= ml+iw) and (mt <= ev.y() <= mt+ih):
			fx, fy = self._display_to_img(ev.x(), ev.y())
			self.canvas_mouse_pos.emit(fx, fy)
		else:
			self.canvas_mouse_out.emit()

	def leaveEvent(self, ev):
		self.canvas_mouse_out.emit()

# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# ROI preview (QPainter-based — see ROICanvas below)
# ---------------------------------------------------------------------------


class ROICanvas(QWidget):
	"""
	ROI preview using the same QPainter pipeline as IRDisplay.
	Always shows the full frame; highlights ROI with dimming + lime border.
	Ticks show true full-frame pixel coordinates.
	"""
	def __init__(self):
		super().__init__()
		self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
		self.setMinimumHeight(140)
		self._img  : Optional[np.ndarray] = None
		self._lut  : Optional[np.ndarray] = None
		self._lo   = 0.0; self._hi = 1.0
		self._roi  = [0,0,0,0]
		self._inv_x = False; self._inv_y = False
		self._MARG_L = 44
		self._MARG_T = 8
		self._MARG_R = 52
		self._MARG_B = 30   # tick margin (pixels)
		self.setMouseTracking(False)

	def update_view(self, img, roi, clim, cmap_name, inv_x, inv_y):
		cmax = _CMAP_MAXCOLS[_CMAP_KEYS.index(cmap_name)] if cmap_name in _CMAP_KEYS else 256
		n    = min(256, cmax)
		lut  = _build_lut(cmap_name, n)
		self._img = img; self._lut = lut
		self._lo, self._hi = clim
		self._roi = list(roi)
		self._inv_x = inv_x; self._inv_y = inv_y
		self.update()

	def paintEvent(self, ev):
		p = QPainter(self)
		p.fillRect(self.rect(), QColor("#e4e4e4"))
		if self._img is None or self._lut is None:
			p.end(); return

		H, W = self._img.shape
		ML = self._MARG_L; MT = self._MARG_T
		MR = self._MARG_R; MB = self._MARG_B
		avail_w = max(self.width()  - ML - MR, 1)
		avail_h = max(self.height() - MT - MB, 1)
		# Preserve aspect ratio — centre the image in available space
		scale   = min(avail_w / W, avail_h / H)
		disp_w  = int(W * scale); disp_h = int(H * scale)
		img_x   = ML + (avail_w - disp_w) // 2
		img_y   = MT + (avail_h - disp_h) // 2

		# Remap to uint8 via LUT
		span = self._hi - self._lo
		if span == 0: span = 1e-9
		normed = np.clip((self._img - self._lo) / span * 255, 0, 255).astype(np.uint8)
		if self._inv_x: normed = normed[:, ::-1]
		if self._inv_y: normed = normed[::-1, :]
		rgb = np.ascontiguousarray(self._lut[normed][:, :, ::-1])
		qi  = QImage(rgb.tobytes(), W, H, W * 3, QImage.Format_RGB888)
		px  = QPixmap.fromImage(qi).scaled(disp_w, disp_h,
										   Qt.IgnoreAspectRatio, Qt.FastTransformation)
		p.drawPixmap(img_x, img_y, px)

		# ROI shading + lime border
		x0, x1, y0, y1 = self._roi
		is_full = (x0 == 0 and y0 == 0 and x1 == W - 1 and y1 == H - 1)
		if not is_full:
			def _ix(col): return img_x + int(col * scale)
			def _iy(row): return img_y + int(row * scale)
			if self._inv_x: _ix = lambda col, _W=W: img_x + int(((_W - 1) - col) * scale)
			if self._inv_y: _iy = lambda row, _H=H: img_y + int(((_H - 1) - row) * scale)
			rx0 = _ix(x0); rx1 = _ix(x1 + 1)
			ry0 = _iy(y0); ry1 = _iy(y1 + 1)
			if rx0 > rx1: rx0, rx1 = rx1, rx0
			if ry0 > ry1: ry0, ry1 = ry1, ry0
			shade = QColor(0, 0, 0, 110)
			for (sx, sy, sw, sh) in [
				(img_x, img_y,        disp_w, ry0 - img_y),
				(img_x, ry1,          disp_w, img_y + disp_h - ry1),
				(img_x, ry0,          rx0 - img_x,             ry1 - ry0),
				(rx1,   ry0,          img_x + disp_w - rx1,    ry1 - ry0),
			]:
				if sw > 0 and sh > 0: p.fillRect(sx, sy, sw, sh, shade)
			p.setPen(QPen(QColor(0, 255, 0), 2)); p.setBrush(Qt.NoBrush)
			p.drawRect(rx0, ry0, rx1 - rx0, ry1 - ry0)

		# Border around image
		p.setPen(QPen(QColor("#888888"), 1)); p.setBrush(Qt.NoBrush)
		p.drawRect(img_x, img_y, disp_w, disp_h)

		# Axis ticks: X on top (or bottom if not inv_y), Y on left (or right if inv_x)
		font = QFont(); font.setPixelSize(11); p.setFont(font)
		p.setPen(QPen(QColor("#1a1a1a"), 1))
		N_TICKS = 5
		for i in range(N_TICKS):
			frac = i / (N_TICKS - 1)
			col  = int(round(frac * (W - 1)))
			row  = int(round(frac * (H - 1)))
			if self._inv_x: col = W - 1 - col
			if self._inv_y: row = H - 1 - row
			wx = img_x + int(frac * disp_w)
			wy = img_y + int(frac * disp_h)
			# X tick
			x_top = self._inv_y
			if x_top:
				p.drawLine(wx, img_y - 5, wx, img_y)
				p.drawText(wx - 20, img_y - 19, 40, 16, Qt.AlignCenter, str(col))
			else:
				p.drawLine(wx, img_y + disp_h, wx, img_y + disp_h + 5)
				p.drawText(wx - 20, img_y + disp_h + 6, 40, 16, Qt.AlignCenter, str(col))
			# Y tick
			y_left = not self._inv_x
			if y_left:
				p.drawLine(img_x - 5, wy, img_x, wy)
				p.drawText(0, wy - 8, img_x - 6, 16, Qt.AlignRight | Qt.AlignVCenter, str(row))
			else:
				p.drawLine(img_x + disp_w, wy, img_x + disp_w + 5, wy)
				p.drawText(img_x + disp_w + 6, wy - 8, MR - 4, 16,
						   Qt.AlignLeft | Qt.AlignVCenter, str(row))
		p.end()


# ---------------------------------------------------------------------------
# Profile window (unchanged from previous version)
# ---------------------------------------------------------------------------

class ProfileWindow(QMainWindow):
	params_changed = pyqtSignal()
	def __init__(self, parent=None):
		super().__init__(parent)
		self.setWindowTitle("Extrapolated Contrast - Temperature Profile")
		self.resize(960,400); self.setVisible(False)
		self._lines=[]; self._ti0=0; self.delta_t=1.0
		cw=QWidget(); self.setCentralWidget(cw)
		vl=QVBoxLayout(cw); vl.setContentsMargins(10,10,10,10); vl.setSpacing(8)
		self.fig_b=Figure(tight_layout=True); self.canvas_b=FigureCanvas(self.fig_b)
		self.ax_b=self.fig_b.add_subplot(111)
		self.ax_b.set_xscale("log"); self.ax_b.set_yscale("log")
		self.ax_b.set_title(r"$\Delta T$ vs. $t$ (bi-log)")
		vl.addWidget(self.canvas_b,stretch=3)
		ctrl=QWidget(); gl=QGridLayout(ctrl); gl.setSpacing(6); vl.addWidget(ctrl)
		def _sl(lo,hi,v,en=True):
			s=QSlider(Qt.Horizontal); s.setRange(lo,hi); s.setValue(v); s.setEnabled(en); return s
		gl.addWidget(QLabel("Init t gross:"),0,0); self.sl_ig=_sl(0,1000,0); gl.addWidget(self.sl_ig,0,1)
		gl.addWidget(QLabel("Init t fine:"), 1,0); self.sl_if=_sl(0,100,0);  gl.addWidget(self.sl_if,1,1)
		self.lbl_it=QLabel("t0=0 s"); gl.addWidget(self.lbl_it,0,2,2,1)
		gl.addWidget(QLabel("Diffusivity gross:"),0,3); self.sl_dg=_sl(0,1000,0,False); gl.addWidget(self.sl_dg,0,4)
		gl.addWidget(QLabel("Diffusivity fine:"), 1,3); self.sl_df=_sl(0,100,10,False); gl.addWidget(self.sl_df,1,4)
		self.lbl_dif=QLabel("a=0.10 mm2/s"); gl.addWidget(self.lbl_dif,0,5,2,1)
		gl.addWidget(QLabel("Plate L:"),2,3); self.sl_L=_sl(0,1000,400,False); gl.addWidget(self.sl_L,2,4)
		self.lbl_L=QLabel("L=4.00 mm"); gl.addWidget(self.lbl_L,2,5)
		gl.addWidget(QLabel("Losses h:"),3,3); self.sl_h=_sl(0,350,100,False); gl.addWidget(self.sl_h,3,4)
		self.lbl_h=QLabel("h=10"); gl.addWidget(self.lbl_h,3,5)
		gl.addWidget(QLabel("Cond. gross:"),2,0); self.sl_cg=_sl(0,2500,0,False); gl.addWidget(self.sl_cg,2,1)
		gl.addWidget(QLabel("Cond. fine:"), 3,0); self.sl_cf=_sl(0,250,20,False); gl.addWidget(self.sl_cf,3,1)
		self.lbl_cond=QLabel("lam=0.20 W/mK"); gl.addWidget(self.lbl_cond,2,2,2,1)
		self.cb_etype=QComboBox(); self.cb_etype.addItems(_EXTRAP_TYPES)
		gl.addWidget(QLabel("Model:"),4,3); gl.addWidget(self.cb_etype,4,4,1,2)
		br=QHBoxLayout()
		def _tb(t,ch=False): b=QPushButton(t); b.setCheckable(True); b.setChecked(ch); return b
		self.btn_clear=QPushButton("Clear"); self.btn_ext=_tb("Ext.")
		self.btn_logeq=_tb("Lin/Log"); self.btn_logsp=_tb("LogSc",True)
		self.btn_norm=_tb("Norm."); self.btn_lock=_tb("Auto",True)
		self.btn_grid=_tb("Grid"); self.btn_leg=_tb("Legend")
		for b in [self.btn_clear,self.btn_ext,self.btn_logeq,self.btn_logsp,
				  self.btn_norm,self.btn_lock,self.btn_grid,self.btn_leg]: br.addWidget(b)
		gl.addLayout(br,5,0,1,6)
		for sl in [self.sl_ig,self.sl_if,self.sl_dg,self.sl_df,
				   self.sl_L,self.sl_h,self.sl_cg,self.sl_cf]:
			sl.valueChanged.connect(self._on_p)
		self.cb_etype.currentIndexChanged.connect(self._on_p)
		self.btn_clear.clicked.connect(self.clear_all)
		self.btn_ext.toggled.connect(self._redraw)
		self.btn_logeq.toggled.connect(lambda on:(self.btn_logeq.setText("LogEq" if on else "LinEq"),self._redraw()))
		self.btn_logsp.toggled.connect(lambda on:(self.ax_b.set_xscale("log" if on else "linear"),
												   self.ax_b.set_yscale("log" if on else "linear"),
												   self.canvas_b.draw_idle()))
		self.btn_lock.toggled.connect(lambda _:self._redraw())
		self.btn_grid.toggled.connect(lambda on:(self.ax_b.grid(on),self.canvas_b.draw_idle()))
		self.btn_leg.toggled.connect(self._redraw)
		self.btn_norm.toggled.connect(lambda _:self.params_changed.emit())

	@property
	def etype(self): return self.cb_etype.currentIndex()+1
	@property
	def L_m(self): return self.sl_L.value()/100*1e-3
	@property
	def diff(self): return (self.sl_dg.value()/10+self.sl_df.value()/100)*1e-6
	@property
	def cond(self): return self.sl_cg.value()/10+self.sl_cf.value()/100
	@property
	def h(self): return self.sl_h.value()/10
	@property
	def t0(self): return (self.sl_ig.value()+self.sl_if.value()/10)/100*self.delta_t
	@property
	def ti0(self): return self._ti0
	@ti0.setter
	def ti0(self,v): self._ti0=int(v); self._redraw()
	def _xdata(self,N): return np.arange(N)*self.delta_t+self.t0+1e-5

	def _on_p(self):
		self.lbl_it.setText(f"t0={self.t0:.4g}s")
		self.lbl_dif.setText(f"a={self.diff*1e6:.3g}mm2/s")
		self.lbl_L.setText(f"L={self.L_m*1e3:.3g}mm")
		self.lbl_h.setText(f"h={self.h:.3g}")
		self.lbl_cond.setText(f"lam={self.cond:.3g}W/mK")
		plate=self.etype>=2; losses=self.etype==3
		for w in [self.sl_dg,self.sl_df,self.sl_L]: w.setEnabled(plate)
		for w in [self.sl_h,self.sl_cg,self.sl_cf]: w.setEnabled(losses)
		self._redraw(); self.params_changed.emit()

	def add_curve(self,xy,T):
		col=['red','blue','darkcyan','purple','darkgreen'][len(self._lines)%5]
		ln,=self.ax_b.plot([],[],'.', color=col,label=f"DT[{xy[0]},{xy[1]}]")
		ex,=self.ax_b.plot([],[],'--',color=col)
		tp,=self.ax_b.plot([],[],'o', color='green')
		fp,=self.ax_b.plot([],[],'x', color='black',ms=12)
		self._lines.append({"xy":xy,"T":T.copy(),"ln":ln,"ex":ex,"tp":tp,"fp":fp})
		self._redraw()

	def clear_all(self):
		for e in self._lines:
			for k in ("ln","ex","tp","fp"):
				try: e[k].remove()
				except: pass
		self._lines.clear(); self.canvas_b.draw_idle()

	def update_frame_marker(self,frame):
		logeq=self.btn_logeq.isChecked()
		for e in self._lines:
			N=len(e["T"]); xd=self._xdata(N); xi=np.full(N,np.nan)
			if frame<N: xi[frame]=xd[frame]
			yi=e["T"].astype(float)
			if logeq:
				xi=np.where(np.isfinite(xi),np.log(np.maximum(xi,1e-30)),np.nan)
				yi=np.log(np.maximum(yi,1e-30))
			e["fp"].set_data(xi,yi)
		self.canvas_b.draw_idle()

	def _redraw(self,_=None):
		logeq=self.btn_logeq.isChecked(); show_ext=self.btn_ext.isChecked()
		for e in self._lines:
			T=e["T"].astype(float); N=len(T); xd=self._xdata(N)
			x=np.log(np.maximum(xd,1e-30)) if logeq else xd
			y=np.log(np.maximum(T,1e-30))  if logeq else T
			e["ln"].set_data(x,y)
			if show_ext and self._ti0<N:
				try:
					ET=extrapolated_temperature(xd[self._ti0],float(T[self._ti0]),xd,
						self.etype,L=self.L_m,diff=self.diff,cond=self.cond,h=self.h)
					ex=np.log(np.maximum(xd,1e-30)) if logeq else xd
					ey=np.log(np.maximum(np.abs(ET),1e-30)) if logeq else np.real(ET)
					e["ex"].set_data(ex,ey)
					xi=np.full(N,np.nan); xi[self._ti0]=x[self._ti0]
					e["tp"].set_data(xi,y)
				except: e["ex"].set_data([],[])
			else: e["ex"].set_data([],[])
		if self._lines:
			if not self.btn_lock.isChecked(): self.ax_b.relim(); self.ax_b.autoscale_view()
			if self.btn_leg.isChecked():
				self.ax_b.legend([e["ln"] for e in self._lines],
								 [e["ln"].get_label() for e in self._lines])
			else:
				lg=self.ax_b.get_legend()
				if lg: lg.remove()
		self.canvas_b.draw_idle()

# ---------------------------------------------------------------------------
# Display panel (left pane, below FramePanel)
# ---------------------------------------------------------------------------

class DisplayPanel(QWidget):
	def __init__(self):
		super().__init__()
		self.setStyleSheet(_STYLE)
		g  = QGroupBox("Display"); gl = QGridLayout(g); gl.setSpacing(7)
		self.chk_inv_x = QCheckBox("Invert X")
		self.chk_inv_y = QCheckBox("Invert Y")
		self.chk_axis  = QCheckBox("Axis"); self.chk_axis.setChecked(True)
		self.chk_grid  = QCheckBox("Grid")
		gl.addWidget(self.chk_inv_x, 0, 0); gl.addWidget(self.chk_inv_y, 0, 1)
		gl.addWidget(self.chk_axis,  1, 0); gl.addWidget(self.chk_grid,  1, 1)
		vl = QVBoxLayout(self); vl.setContentsMargins(6, 4, 6, 4)
		vl.addWidget(g)

# ---------------------------------------------------------------------------
# Frame panel (left pane, between ROI preview and action bar)
# ---------------------------------------------------------------------------

class FramePanel(QWidget):
	def __init__(self, n_frames=0):
		super().__init__()
		self.setStyleSheet(_STYLE)
		g = QGroupBox("Frame"); gl = QGridLayout(g); gl.setSpacing(7)
		self.lbl_frame = QLabel(f"1 / {max(n_frames,1)}")
		self.lbl_frame.setAlignment(Qt.AlignCenter)
		gl.addWidget(self.lbl_frame,0,0,1,4)
		self.sl_frame = QSlider(Qt.Horizontal)
		self.sl_frame.setRange(1,max(n_frames,1)); self.sl_frame.setValue(1)
		self.sl_frame.setEnabled(n_frames>1)
		gl.addWidget(self.sl_frame,1,0,1,4)
		self.btn_first = QPushButton("⏮"); self.btn_first.setFixedSize(38,30)
		self.btn_bk1   = QPushButton("◀");  self.btn_bk1.setFixedSize(38,30)
		self.btn_fw1   = QPushButton("▶");  self.btn_fw1.setFixedSize(38,30)
		self.btn_last  = QPushButton("⏭");  self.btn_last.setFixedSize(38,30)
		self.btn_bk1.setAutoRepeat(True); self.btn_bk1.setAutoRepeatDelay(400); self.btn_bk1.setAutoRepeatInterval(80)
		self.btn_fw1.setAutoRepeat(True); self.btn_fw1.setAutoRepeatDelay(400); self.btn_fw1.setAutoRepeatInterval(80)
		nr = QHBoxLayout(); nr.setSpacing(4)
		for b in [self.btn_first,self.btn_bk1,self.btn_fw1,self.btn_last]: nr.addWidget(b)
		gl.addLayout(nr,2,0,1,4)
		vl = QVBoxLayout(self); vl.setContentsMargins(6,6,6,6)
		vl.addWidget(g)

	def set_enabled(self, on): self.setEnabled(on)
	def set_frame_label(self, f1, total): self.lbl_frame.setText(f"{f1} / {total}")

	@property
	def current_frame(self): return self.sl_frame.value()-1

# ---------------------------------------------------------------------------
# Control panel (right pane — no Frame section)
# ---------------------------------------------------------------------------

class ControlPanel(QWidget):
	def __init__(self):
		super().__init__()
		self.setStyleSheet(_STYLE); self.setMinimumWidth(256)
		root = QVBoxLayout(self); root.setContentsMargins(8,8,8,8); root.setSpacing(10)

		# Spatial Filter
		g3 = QGroupBox("Spatial Filter"); gl3 = QGridLayout(g3); gl3.setSpacing(7)
		self.btn_noflt  = QPushButton("None");     self.btn_noflt.setCheckable(True);  self.btn_noflt.setChecked(True)
		self.btn_gauss  = QPushButton("Gaussian"); self.btn_gauss.setCheckable(True)
		self.btn_median = QPushButton("Median");   self.btn_median.setCheckable(True)
		self._flt_grp = QButtonGroup(); self._flt_grp.setExclusive(True)
		for b in [self.btn_noflt,self.btn_gauss,self.btn_median]: self._flt_grp.addButton(b)
		gl3.addWidget(self.btn_noflt,0,0); gl3.addWidget(self.btn_gauss,0,1); gl3.addWidget(self.btn_median,0,2)
		self.lbl_var = QLabel("Variance: 0.85"); self.lbl_var.setEnabled(False)
		self.sl_var  = QSlider(Qt.Horizontal); self.sl_var.setRange(5,500); self.sl_var.setValue(85); self.sl_var.setEnabled(False)
		gl3.addWidget(self.lbl_var,1,0,1,3); gl3.addWidget(self.sl_var,2,0,1,3)
		root.addWidget(g3)

		# Colour / Clamping
		g4 = QGroupBox("Colour / Clamping"); gl4 = QGridLayout(g4); gl4.setSpacing(7)
		self.cb_cmap = QComboBox(); self.cb_cmap.addItems(_CMAP_LABELS)
		self.cb_cmap.setCurrentIndex(_CMAP_KEYS.index("jet"))
		gl4.addWidget(self.cb_cmap,0,0,1,2)
		self.lbl_ncolors = QLabel("Colors: 256 / 256")
		gl4.addWidget(self.lbl_ncolors,1,0,1,2)
		self.sl_ncolors = QSlider(Qt.Horizontal); self.sl_ncolors.setRange(2,256); self.sl_ncolors.setValue(256)
		gl4.addWidget(self.sl_ncolors,2,0,1,2)
		# 2x2 clamp buttons
		self.btn_clamp_img = QPushButton("Fit Frame"); self.btn_clamp_img.setCheckable(True); self.btn_clamp_img.setChecked(True)
		self.btn_clamp_seq = QPushButton("Fit Seq");   self.btn_clamp_seq.setCheckable(True)
		self.btn_clamp_roi = QPushButton("Fit ROI");   self.btn_clamp_roi.setCheckable(True)
		self.btn_clamp_man = QPushButton("Manual");    self.btn_clamp_man.setCheckable(True)
		self._clamp_grp = QButtonGroup(); self._clamp_grp.setExclusive(True)
		for b in [self.btn_clamp_img,self.btn_clamp_seq,self.btn_clamp_roi,self.btn_clamp_man]:
			self._clamp_grp.addButton(b)
		gl4.addWidget(self.btn_clamp_img,3,0); gl4.addWidget(self.btn_clamp_seq,3,1)
		gl4.addWidget(self.btn_clamp_roi,4,0); gl4.addWidget(self.btn_clamp_man,4,1)
		gl4.addWidget(QLabel("High"),5,0)
		self.edt_hi = QLineEdit("0.0"); self.edt_hi.setValidator(QDoubleValidator())
		gl4.addWidget(self.edt_hi,5,1)
		self.sl_hi = QSlider(Qt.Horizontal); gl4.addWidget(self.sl_hi,6,0,1,2)
		gl4.addWidget(QLabel("Low"),7,0)
		self.edt_lo = QLineEdit("0.0"); self.edt_lo.setValidator(QDoubleValidator())
		gl4.addWidget(self.edt_lo,7,1)
		self.sl_lo = QSlider(Qt.Horizontal); gl4.addWidget(self.sl_lo,8,0,1,2)
		self.chk_lock = QCheckBox("Lock all"); gl4.addWidget(self.chk_lock,9,0,1,2)
		self._manual_widgets = [self.edt_hi,self.sl_hi,self.edt_lo,self.sl_lo]
		for w in self._manual_widgets: w.setEnabled(False)
		root.addWidget(g4)

		# Processing
		g5 = QGroupBox("Processing"); gl5 = QGridLayout(g5); gl5.setSpacing(7)
		def _p(t):
			b = QPushButton(t); b.setCheckable(True); return b

		self.btn_temp    = _p("Temperature"); self.btn_temp.setChecked(True)
		self.btn_ftamp   = _p("FT Amplitude")
		self.btn_ftphase = _p("FT Phase")
		self.btn_extrap  = _p("Extrap. CT")
		self.btn_corr    = _p("Correlation")
		self.btn_setupct = _p("Setup CT")
		# New algorithms
		self.btn_ca      = _p("Abs. Contrast")
		self.btn_dac     = _p("Diff. AC")
		self.btn_rx      = _p("RX Detector")
		self.btn_skq     = _p("SKQ")
		self.btn_pca_m   = _p("PCA M-Method")
		self.btn_tsr     = _p("TSR")
		self.btn_wavelet = _p("Wavelet")

		self._proc_grp = QButtonGroup(); self._proc_grp.setExclusive(True)
		for b in [self.btn_temp, self.btn_ftamp, self.btn_ftphase,
				  self.btn_extrap, self.btn_corr,
				  self.btn_ca, self.btn_dac, self.btn_rx, self.btn_skq,
				  self.btn_pca_m, self.btn_tsr, self.btn_wavelet]:
			self._proc_grp.addButton(b)

		# SKQ sub-mode selector (shown only when SKQ active)
		self.cb_skq = QComboBox()
		self.cb_skq.addItems(["Kurtosis", "Skewness", "5th Moment"])
		self.cb_skq.setEnabled(False)

		r = 0
		gl5.addWidget(self.btn_temp,    r, 0, 1, 2); r += 1
		gl5.addWidget(self.btn_ftamp,   r, 0); gl5.addWidget(self.btn_ftphase, r, 1); r += 1
		gl5.addWidget(self.btn_extrap,  r, 0, 1, 2); r += 1
		gl5.addWidget(self.btn_corr,    r, 0); gl5.addWidget(self.btn_setupct,  r, 1); r += 1
		gl5.addWidget(self.btn_ca,      r, 0); gl5.addWidget(self.btn_dac,     r, 1); r += 1
		gl5.addWidget(self.btn_rx,      r, 0); gl5.addWidget(self.btn_skq,     r, 1); r += 1
		gl5.addWidget(self.cb_skq,      r, 0, 1, 2); r += 1
		gl5.addWidget(self.btn_pca_m,   r, 0); gl5.addWidget(self.btn_tsr,     r, 1); r += 1
		gl5.addWidget(self.btn_wavelet, r, 0, 1, 2); r += 1

		sep = QFrame(); sep.setFrameShape(QFrame.HLine); gl5.addWidget(sep, r, 0, 1, 2); r += 1
		self.btn_set_roi   = QPushButton("Set ROI")
		self.btn_reset_roi = QPushButton("Reset ROI")
		gl5.addWidget(self.btn_set_roi,   r, 0)
		gl5.addWidget(self.btn_reset_roi, r, 1)
		root.addWidget(g5)

		self._data_groups = [g3, g4, g5]
		root.addStretch(1)

	@property
	def variance(self): v=self.sl_var.value()/100.; return round(v*100/5)/100*5
	@property
	def cmap_name(self): return _CMAP_KEYS[self.cb_cmap.currentIndex()]
	@property
	def n_colors(self): return self.sl_ncolors.value()
	@property
	def clamp_mode(self):
		if self.btn_clamp_img.isChecked(): return CLAMP_IMG
		if self.btn_clamp_seq.isChecked(): return CLAMP_SEQ
		if self.btn_clamp_roi.isChecked(): return CLAMP_ROI
		return CLAMP_MAN

	def set_data_enabled(self, on):
		for w in self._data_groups: w.setEnabled(on)
	def set_manual_enabled(self, on):
		for w in self._manual_widgets: w.setEnabled(on)
	def set_lock(self, on):
		locked = [self.cb_cmap,self.sl_ncolors,
				  self.btn_clamp_img,self.btn_clamp_seq,
				  self.btn_clamp_roi,self.btn_clamp_man]+self._manual_widgets
		for w in locked: w.setEnabled(not on)
	def set_clamp(self, lo, hi):
		self.edt_lo.setText(f"{lo:.4g}"); self.edt_hi.setText(f"{hi:.4g}")
		self._cb_lo=lo; self._cb_hi=hi
		for sl in [self.sl_lo,self.sl_hi]: sl.setRange(0,10000)
		self.sl_lo.setValue(0); self.sl_hi.setValue(10000)
	def get_clamp_sliders(self):
		lo=getattr(self,'_cb_lo',0.); hi=getattr(self,'_cb_hi',1.); span=hi-lo
		return lo+self.sl_lo.value()/10000.*span, lo+self.sl_hi.value()/10000.*span
	def get_clamp_edits(self):
		try: return float(self.edt_lo.text()), float(self.edt_hi.text())
		except ValueError: return None,None

# ---------------------------------------------------------------------------
# Action bar (left pane bottom - 2 column grid, fixed height)
# ---------------------------------------------------------------------------

class ActionBar(QWidget):
	def __init__(self):
		super().__init__()
		self.setStyleSheet(_STYLE)
		self.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Fixed)
		gl = QGridLayout(self); gl.setContentsMargins(8,6,8,10); gl.setSpacing(6)
		self.btn_open    = QPushButton("Open file")
		self.btn_stream  = QPushButton("Open stream")
		self.btn_save    = QPushButton("Save session")
		self.btn_png     = QPushButton("Export PNG")
		self.btn_capture = QPushButton("Capture sequence")
		self.btn_exit    = QPushButton("Exit")
		for b in [self.btn_open, self.btn_stream, self.btn_save,
				  self.btn_png, self.btn_capture, self.btn_exit]:
			b.setMinimumHeight(28)
		btns = [self.btn_open, self.btn_stream, self.btn_save,
				self.btn_png, self.btn_capture, self.btn_exit]
		for i, b in enumerate(btns):
			gl.addWidget(b, i // 2, i % 2)
		self._data_btns = [self.btn_save, self.btn_png, self.btn_capture]
		self.set_data_enabled(False)

	def set_data_enabled(self, on):
		for b in self._data_btns: b.setEnabled(on)
	def set_file_open(self, on):
		self.btn_open.setText("Close file" if on else "Open file")
	def set_streaming(self, on):
		self.btn_stream.setText("Close stream" if on else "Open stream")

# ---------------------------------------------------------------------------
# Capture dialog
# ---------------------------------------------------------------------------

class CaptureDialog(QDialog):
	def __init__(self, parent=None):
		super().__init__(parent)
		self.setWindowTitle("Capture sequence")
		self.setStyleSheet(_STYLE)
		self.setMinimumWidth(320)
		vl = QVBoxLayout(self); vl.setContentsMargins(16,16,16,12); vl.setSpacing(10)

		vl.addWidget(QLabel("Number of samples:"))
		self.spin_samples = QSpinBox()
		self.spin_samples.setRange(1,10000); self.spin_samples.setValue(50)
		vl.addWidget(self.spin_samples)

		vl.addWidget(QLabel("Total capture time (seconds):"))
		self.spin_time = QDoubleSpinBox()
		self.spin_time.setRange(0.1,3600.0); self.spin_time.setValue(10.0)
		self.spin_time.setDecimals(1); self.spin_time.setSingleStep(0.5)
		vl.addWidget(self.spin_time)

		self.lbl_rate = QLabel("")
		vl.addWidget(self.lbl_rate)

		bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
		bb.accepted.connect(self._on_accept); bb.rejected.connect(self.reject)
		vl.addWidget(bb)

		self.spin_samples.valueChanged.connect(self._update_rate)
		self.spin_time.valueChanged.connect(self._update_rate)
		self._update_rate()

	def _update_rate(self):
		s = self.spin_samples.value(); t = self.spin_time.value()
		rate = s/t if t>0 else 0
		warning = "  *** exceeds 25 fps camera limit!" if rate>25 else ""
		self.lbl_rate.setText(f"Rate: {rate:.2f} fps{warning}")
		self.lbl_rate.setStyleSheet("color:red;" if rate>25 else "color:#1a1a1a;")

	def _on_accept(self):
		s = self.spin_samples.value(); t = self.spin_time.value()
		if s/t > 25.01:
			QMessageBox.warning(self,"Rate too high",
				"samples/time cannot exceed 25 fps (camera limit).\n"
				"Reduce samples or increase time.")
			return
		self.accept()

	@property
	def n_samples(self): return self.spin_samples.value()
	@property
	def total_time(self): return self.spin_time.value()

# ---------------------------------------------------------------------------
# Main application
# ---------------------------------------------------------------------------

class IrView:
	def __init__(self):
		self.app = QApplication.instance() or QApplication(sys.argv)
		self.seq  : Optional[np.ndarray] = None
		self.H=self.W=self.N=0
		self.roi  = [0,0,0,0]
		self._fft_cache  = self._fft_type = None
		self._corr_cache = None
		self._pct_U      = self._pct_S = None  # kept for session compat
		self._ca_cache   : Optional[np.ndarray] = None
		self._dac_cache  : Optional[np.ndarray] = None
		self._rx_cache   : Optional[np.ndarray] = None
		self._skq_cache  : Optional[dict]       = None
		self._skq_key    : str                  = "kurtosis"
		self._pca_cache  : Optional[np.ndarray] = None
		self._pca_S      : Optional[np.ndarray] = None
		self._tsr_synth  : Optional[np.ndarray] = None
		self._haar_ca    : Optional[np.ndarray] = None
		self._haar_cd    : Optional[np.ndarray] = None
		self._cold_frame : Optional[np.ndarray] = None
		self._ca_ref_px  : tuple = (0, 0)
		self._dac_tprime : int   = 1
		self._tsr_degree : int   = 5
		self._roi_selecting = False
		self._roi_pt1 : Optional[tuple] = None
		self._roi_mpl_cid = self._roi_move_cid = None
		self._loading = False
		self._source_path = ""
		self._mode = MODE_IDLE
		self._cam = CameraReceiver()
		self._cam_timer : Optional[QTimer] = None
		self._cam_fidx_last = -1
		self._cap_buf   : list = []
		self._cap_timer : Optional[QTimer] = None
		self._cap_n     = 0
		self._cap_path  = ""
		self._hover_x   = None
		self._hover_y   = None
		self._hover_z   = 0.0
		self._build_ui()

	# -- layout ---------------------------------------------------------------

	def _build_ui(self):
		self.main_win = QMainWindow()
		self.main_win.setWindowTitle("IR View")

		self._sb = self.main_win.statusBar()
		self._sb.setStyleSheet(
			"QStatusBar{background:#dcdcdc;color:#111111;"
			"font-size:14px;padding:3px 10px;border-top:1px solid #b0b0b0;}")
		self._sb.showMessage("No file loaded.")

		# Left column: ROI preview | frame panel | display panel | action bar
		self.roi_canvas    = ROICanvas()
		self.frame_panel   = FramePanel(0)
		self.display_panel = DisplayPanel()
		self.action_bar    = ActionBar()

		left_col = QWidget(); left_col.setStyleSheet(_STYLE)
		left_vl  = QVBoxLayout(left_col)
		left_vl.setContentsMargins(0, 0, 0, 0); left_vl.setSpacing(0)
		left_vl.addWidget(self.roi_canvas, stretch=1)
		left_vl.addWidget(self.frame_panel, stretch=0)
		left_vl.addWidget(self.display_panel, stretch=0)
		sep1 = QFrame(); sep1.setFrameShape(QFrame.HLine)
		left_vl.addSpacing(6)
		left_vl.addWidget(self.action_bar, stretch=0)

		# Centre: main display
		self.canvas = IRDisplay()
		canvas_wrap = QWidget()
		ch = QHBoxLayout(canvas_wrap); ch.setContentsMargins(6,0,0,0)
		ch.addWidget(self.canvas)

		# Right: control panel
		self.ctrl = ControlPanel()
		self._ctrl_scroll = QScrollArea()
		self._ctrl_scroll.setWidget(self.ctrl); self._ctrl_scroll.setWidgetResizable(True)
		self._ctrl_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
		self._ctrl_scroll.setMinimumWidth(272)

		h_split = QSplitter(Qt.Horizontal)
		h_split.addWidget(left_col); h_split.addWidget(canvas_wrap); h_split.addWidget(self._ctrl_scroll)
		h_split.setStretchFactor(0,1); h_split.setStretchFactor(1,4); h_split.setStretchFactor(2,1)
		h_split.setSizes([250,830,280])
		self._h_split = h_split
		self.main_win.setCentralWidget(h_split)

		self.profile_win = ProfileWindow(self.main_win)

		QShortcut(QKeySequence(Qt.Key_Escape), self.main_win).activated.connect(self._cancel_roi)

		self._connect_static_signals()
		self._connect_ctrl_signals()
		self.ctrl.set_data_enabled(False)
		self.display_panel.setEnabled(False)
		self.frame_panel.set_enabled(False)
		self.main_win.resize(1360, 800)

	def _connect_static_signals(self):
		ab = self.action_bar
		ab.btn_open.clicked.connect(self._open_dialog)
		ab.btn_stream.clicked.connect(self._toggle_stream)
		ab.btn_save.clicked.connect(self._save_session)
		ab.btn_png.clicked.connect(self._export_png)
		ab.btn_capture.clicked.connect(self._start_capture_dialog)
		ab.btn_exit.clicked.connect(self.main_win.close)
		self.canvas.canvas_mouse_pos.connect(self._on_mouse_pos)
		self.canvas.canvas_mouse_out.connect(self._on_mouse_out)
		self.profile_win.params_changed.connect(self._refresh)

	def _connect_frame_signals(self):
		fp = self.frame_panel
		fp.sl_frame.valueChanged.connect(self._on_frame)
		fp.btn_first.clicked.connect(lambda: fp.sl_frame.setValue(1))
		fp.btn_last.clicked.connect( lambda: fp.sl_frame.setValue(self.N))
		fp.btn_bk1.clicked.connect(  lambda: fp.sl_frame.setValue(max(1,fp.sl_frame.value()-1)))
		fp.btn_fw1.clicked.connect(  lambda: fp.sl_frame.setValue(min(self.N,fp.sl_frame.value()+1)))

	def _connect_ctrl_signals(self):
		c = self.ctrl
		d = self.display_panel
		d.chk_axis.toggled.connect(self._on_axis)
		d.chk_grid.toggled.connect(self.canvas.toggle_grid)
		d.chk_inv_x.toggled.connect(lambda on: (self.canvas.set_invert_x(on), self._push_roi_preview()))
		d.chk_inv_y.toggled.connect(lambda on: (self.canvas.set_invert_y(on), self._push_roi_preview()))
		c.btn_noflt.toggled.connect(self._on_filter)
		c.btn_gauss.toggled.connect(self._on_filter)
		c.btn_median.toggled.connect(self._on_filter)
		c.sl_var.valueChanged.connect(self._on_variance)
		c.cb_cmap.currentIndexChanged.connect(self._on_cmap)
		c.sl_ncolors.valueChanged.connect(self._on_ncolors)
		for b in [c.btn_clamp_img,c.btn_clamp_seq,c.btn_clamp_roi,c.btn_clamp_man]:
			b.toggled.connect(self._on_clamp_mode)
		c.sl_hi.valueChanged.connect(self._on_clamp_slider)
		c.sl_lo.valueChanged.connect(self._on_clamp_slider)
		c.edt_hi.returnPressed.connect(self._apply_clamp_edits)
		c.edt_lo.returnPressed.connect(self._apply_clamp_edits)
		c.chk_lock.toggled.connect(c.set_lock)
		c._proc_grp.buttonToggled.connect(
			lambda btn, on: (on and not self._loading) and self._on_proc_changed(btn))
		c.cb_skq.currentIndexChanged.connect(
			lambda _: self._skq_cache and self._refresh())
		c.btn_setupct.toggled.connect(self._toggle_setupct)
		c.btn_set_roi.clicked.connect(self._start_roi)
		c.btn_reset_roi.clicked.connect(self._reset_roi)
		try: self.canvas.canvas_clicked.disconnect()
		except: pass
		self.canvas.canvas_clicked.connect(self._on_canvas_click)

	# -- processing guard -----------------------------------------------------

	def _on_proc_changed(self, btn):
		c = self.ctrl
		multi = [c.btn_ftamp, c.btn_ftphase, c.btn_corr,
				 c.btn_extrap, c.btn_ca, c.btn_dac, c.btn_rx, c.btn_skq,
				 c.btn_pca_m, c.btn_tsr, c.btn_wavelet]
		if btn in multi and self.N <= 1:
			c.btn_temp.setChecked(True)
			QMessageBox.information(self.main_win, "Single frame",
				"This mode requires more than one frame.")
			return
		# SKQ sub-mode selector enabled only when SKQ active
		c.cb_skq.setEnabled(btn is c.btn_skq)
		# Clear relevant caches on mode switch
		self._corr_cache = None
		self._ca_cache   = self._dac_cache = None
		self._rx_cache   = self._skq_cache = None
		self._pca_cache  = self._pca_S    = None
		self._tsr_synth  = None
		self._haar_ca    = self._haar_cd  = None
		# CA/DAC need a cold frame — use frame 0 if not set
		if btn in (c.btn_ca, c.btn_dac) and self._cold_frame is None:
			self._cold_frame = self.seq[:, :, 0].copy()
			self._sb.showMessage("Cold frame set to frame 0.")
		# Force Fit Frame clamping on mode change so range adapts
		c.btn_clamp_img.setChecked(True)
		c.set_manual_enabled(False)
		# Frame navigation irrelevant for static (single-image) results
		static_modes = [c.btn_rx, c.btn_skq, c.btn_corr]
		self.frame_panel.set_enabled(btn not in static_modes and self.N > 1)
		self._refresh()

	# -- load -----------------------------------------------------------------

	def _open_dialog(self):
		if self._mode == MODE_FILE:
			# Close current file, return to idle
			self.seq = None; self.H = self.W = self.N = 0
			self._mode = MODE_IDLE; self._source_path = ""
			self.action_bar.set_file_open(False)
			self.action_bar.set_data_enabled(False)
			self.frame_panel.set_enabled(False)
			self.ctrl.set_data_enabled(False)
			self.canvas._img = None; self.canvas.update()
			self.main_win.setWindowTitle("IR View")
			self._sb.showMessage("No file loaded.")
			return
		path, _ = QFileDialog.getOpenFileName(
			self.main_win, "Open file", "", "MATLAB & sessions (*.mat *.irv)")
		if not path: return
		if self._mode == MODE_STREAM: self._close_stream()
		if path.endswith(".irv"): self._load_session(path)
		else: self._load_mat(path)

	def _load_mat(self, path):
		try: seq = load_mat(path)
		except Exception as e:
			QMessageBox.critical(self.main_win,"Load error",str(e)); return
		self._source_path = path; self._mode = MODE_FILE
		self._ingest(seq)
		self.action_bar.set_file_open(True)
		self.main_win.setWindowTitle(f"IR View - {os.path.basename(path)}")

	def _ingest(self, seq):
		self.seq = _ensure3d(seq.astype(np.float32))
		self.H,self.W,self.N = self.seq.shape
		self.roi = [0,self.W-1,0,self.H-1]
		self._fft_cache  = self._fft_type = None
		self._corr_cache = None
		self._ca_cache   = self._dac_cache = None
		self._rx_cache   = self._skq_cache = None
		self._pca_cache  = self._pca_S    = None
		self._tsr_synth  = None
		self._haar_ca    = self._haar_cd  = None
		self._cold_frame = None
		self._roi_selecting = False; self._roi_pt1 = None
		self.canvas.roi_clear_overlay()

		# Rebuild frame panel (remove old, insert new)
		old_fp = self.frame_panel
		self.frame_panel = FramePanel(self.N)
		left_col = self._h_split.widget(0)
		lv = left_col.layout()
		lv.removeWidget(old_fp); old_fp.setParent(None)
		# Insert after roi_canvas (index 0), before display_panel
		lv.insertWidget(1, self.frame_panel, stretch=0)
		self._connect_frame_signals()

		self.frame_panel.set_enabled(self.N>1)
		self.frame_panel.set_frame_label(1,self.N)
		self.ctrl.set_data_enabled(True)
		self.action_bar.set_data_enabled(True)
		self.display_panel.setEnabled(True)
		self.profile_win.delta_t = 1.0
		self._refresh()
		self._do_fit_img()
		self._update_idle_status()

	def _update_idle_status(self):
		if self._mode == MODE_STREAM:
			fn = "Live stream"
		elif self._source_path:
			fn = os.path.basename(self._source_path)
		else:
			fn = ""
		parts = []
		if fn: parts.append(fn)
		if self.W and self.H:
			parts.append(f"{self.W} x {self.H}")
			parts.append(f"{self.N} frame{'s' if self.N != 1 else ''}")
		if self._hover_x is not None:
			parts.append(f"x={self._hover_x}   y={self._hover_y}   z={self._hover_z:.4g}")
		self._sb.showMessage("    |    ".join(parts) if parts else "No file loaded.")

	# -- stream ---------------------------------------------------------------

	def _toggle_stream(self):
		if self._mode == MODE_STREAM:
			self._close_stream()
		else:
			self._open_stream_dialog()

	def _open_stream_dialog(self):
		from PyQt5.QtWidgets import QInputDialog
		host = os.environ.get("IRVIEW_CAM_HOST","127.0.0.1")
		port = os.environ.get("IRVIEW_CAM_PORT","54321")
		txt,ok = QInputDialog.getText(self.main_win,"Camera stream",
									  "Host:port:",text=f"{host}:{port}")
		if not ok or not txt.strip(): return
		try:
			hs,ps = txt.strip().rsplit(":",1); pi = int(ps)
		except Exception:
			QMessageBox.warning(self.main_win,"Bad address","Expected host:port"); return
		if self._mode == MODE_STREAM: self._close_stream()
		self._cam.start(hs,pi)
		QTimer.singleShot(600, self._check_cam)

	def _check_cam(self):
		if not self._cam.is_running():
			QMessageBox.warning(self.main_win,"Stream",
								"Could not connect: "+self._cam.error); return
		self._mode = MODE_STREAM
		self._source_path = ""
		self.action_bar.set_file_open(False)
		self.action_bar.set_streaming(True)
		self._cam_timer = QTimer(); self._cam_timer.setInterval(33)
		self._cam_timer.timeout.connect(self._poll_cam)
		self._cam_timer.start()
		self.main_win.setWindowTitle("IR View - Live stream")
		self._sb.showMessage("Live stream active.")

	def _poll_cam(self):
		if not self._cam.is_running():
			self._close_stream(); return
		frame,fidx,ts = self._cam.get_latest()
		if frame is None or fidx == self._cam_fidx_last: return
		self._cam_fidx_last = fidx
		H,W = frame.shape
		if H!=self.H or W!=self.W or self.seq is None:
			self.seq = frame[:,:,np.newaxis].copy()
			self.H,self.W,self.N = H,W,1
			self.roi = [0,W-1,0,H-1]
			self._fft_cache=self._fft_type=None
			self._corr_cache=self._pct_U=self._pct_S=None
			self.ctrl.set_data_enabled(True)
			self.display_panel.setEnabled(True)
			self.frame_panel.set_enabled(False)
			self.action_bar.set_data_enabled(True)
			for b in [self.ctrl.btn_ftamp, self.ctrl.btn_ftphase,
					  self.ctrl.btn_corr,
					  self.ctrl.btn_extrap, self.ctrl.btn_ca,
					  self.ctrl.btn_dac,   self.ctrl.btn_rx,
					  self.ctrl.btn_skq,   self.ctrl.btn_pca_m,
					  self.ctrl.btn_tsr,   self.ctrl.btn_wavelet]:
				b.setEnabled(False)
		else:
			self.seq[:,:,0] = frame
		self._refresh()
		self._update_idle_status()

	def _close_stream(self):
		if self._cam_timer: self._cam_timer.stop(); self._cam_timer = None
		self._cam.stop(); self._mode = MODE_IDLE
		self.action_bar.set_streaming(False)
		self.main_win.setWindowTitle("IR View")
		self._sb.showMessage("Stream closed.")

	# -- frame pipeline -------------------------------------------------------

	def _current_idx(self):
		return self.frame_panel.current_frame

	def _compute_frame(self, idx, crop=True):
		if self.seq is None: return None
		idx = int(np.clip(idx, 0, self.N - 1))
		c = self.ctrl; r = self.roi

		# ---- FT amplitude / phase -------------------------------------------
		if c.btn_ftamp.isChecked() or c.btn_ftphase.isChecked():
			arr = self._get_fft(1 if c.btn_ftamp.isChecked() else 2)
			if arr is None: return None
			img = arr[:, :, idx].astype(float)

		# ---- Extrapolated contrast -------------------------------------------
		elif c.btn_extrap.isChecked():
			pw  = self.profile_win; eps = 1e-5
			t   = np.arange(self.N) * pw.delta_t + pw.t0 + eps
			sf  = self.seq.astype(float)
			try:
				ET  = extrapolated_temperature(
					float(t[pw.ti0]), sf[:, :, pw.ti0], float(t[idx]),
					pw.etype, L=pw.L_m, diff=pw.diff, cond=pw.cond, h=pw.h)
				img = np.real(sf[:, :, idx] - ET)
				if pw.btn_norm.isChecked():
					ref = sf[:, :, pw.ti0]; ref[ref == 0] = np.nan; img = img / ref
			except Exception as ex:
				self._sb.showMessage(f"Extrap: {ex}")
				img = self.seq[:, :, idx].astype(float)

		# ---- Correlated contrast --------------------------------------------
		elif c.btn_corr.isChecked():
			if self._corr_cache is None:
				self._sb.showMessage("Computing Correlation...")
				QApplication.processEvents()
				res = correlated_contrast(self.seq)
				if res is None:
					QMessageBox.warning(self.main_win, "Memory", "Not enough memory.")
					c.btn_temp.setChecked(True); return None
				self._corr_cache = res
				self._sb.showMessage("")
			img = self._corr_cache.astype(float)

		# ---- Absolute Contrast ----------------------------------------------
		elif c.btn_ca.isChecked():
			if self._ca_cache is None:
				cold = self._cold_frame if self._cold_frame is not None \
					   else self.seq[:, :, 0]
				self._sb.showMessage("Computing Absolute Contrast...")
				QApplication.processEvents()
				try:
					self._ca_cache = absolute_contrast(
						self.seq, cold, self._ca_ref_px[0], self._ca_ref_px[1])
				except Exception as ex:
					QMessageBox.warning(self.main_win, "CA error", str(ex))
					c.btn_temp.setChecked(True); return None
				finally:
					self._sb.showMessage("")
			img = self._ca_cache[:, :, idx].astype(float)

		# ---- Differential Absolute Contrast ---------------------------------
		elif c.btn_dac.isChecked():
			if self._dac_cache is None:
				cold = self._cold_frame if self._cold_frame is not None \
					   else self.seq[:, :, 0]
				tp   = min(self._dac_tprime, self.N - 1)
				self._sb.showMessage("Computing DAC...")
				QApplication.processEvents()
				try:
					self._dac_cache = differential_absolute_contrast(
						self.seq, cold, tp)
				except Exception as ex:
					QMessageBox.warning(self.main_win, "DAC error", str(ex))
					c.btn_temp.setChecked(True); return None
				finally:
					self._sb.showMessage("")
			dac_N = self._dac_cache.shape[2]
			img   = self._dac_cache[:, :, min(idx, dac_N - 1)].astype(float)

		# ---- RX detector ----------------------------------------------------
		elif c.btn_rx.isChecked():
			if self._rx_cache is None:
				self._sb.showMessage("Computing RX detector...")
				QApplication.processEvents()
				try:
					self._rx_cache = rx_detector(self.seq)
				except Exception as ex:
					QMessageBox.warning(self.main_win, "RX error", str(ex))
					c.btn_temp.setChecked(True); return None
				finally:
					self._sb.showMessage("")
			img = self._rx_cache.astype(float)

		# ---- SKQ (kurtosis / skewness / RX / 5th moment) -------------------
		elif c.btn_skq.isChecked():
			if self._skq_cache is None:
				self._sb.showMessage("Computing SKQ statistics...")
				QApplication.processEvents()
				try:
					self._skq_cache = skq_stats(self.seq)
				except Exception as ex:
					QMessageBox.warning(self.main_win, "SKQ error", str(ex))
					c.btn_temp.setChecked(True); return None
				finally:
					self._sb.showMessage("")
			key_map = {0: "kurtosis", 1: "skewness", 2: "rx", 3: "quinto"}
			self._skq_key = key_map.get(c.cb_skq.currentIndex(), "kurtosis")
			img = self._skq_cache[self._skq_key].astype(float)

		# ---- PCA M-Method ---------------------------------------------------
		elif c.btn_pca_m.isChecked():
			if self._pca_cache is None:
				self._sb.showMessage("Computing PCA M-Method...")
				QApplication.processEvents()
				try:
					self._pca_cache, self._pca_S, _ = pca_mmethod(self.seq)
				except Exception as ex:
					QMessageBox.warning(self.main_win, "PCA error", str(ex))
					c.btn_temp.setChecked(True); return None
				finally:
					self._sb.showMessage("")
			nc  = self._pca_cache.shape[2]
			k   = min(idx, nc - 1)
			img = self._pca_cache[:, :, k].astype(float)
			wt  = 100 * self._pca_S[k] if k < len(self._pca_S) else 0
			c.btn_pca_m.setText(f"PCA {wt:.1f}%")

		# ---- TSR (synthetic reconstruction) ---------------------------------
		elif c.btn_tsr.isChecked():
			if self._tsr_synth is None:
				cold = self._cold_frame if self._cold_frame is not None \
					   else self.seq[:, :, 0]
				self._sb.showMessage("Computing TSR...")
				QApplication.processEvents()
				try:
					_, self._tsr_synth = tsr_polyfit(
						self.seq, cold, self._tsr_degree)
				except Exception as ex:
					QMessageBox.warning(self.main_win, "TSR error", str(ex))
					c.btn_temp.setChecked(True); return None
				finally:
					self._sb.showMessage("")
			img = self._tsr_synth[:, :, idx].astype(float)

		# ---- Haar wavelet (approximation coefficients) ----------------------
		elif c.btn_wavelet.isChecked():
			if self._haar_ca is None:
				self._sb.showMessage("Computing Haar DWT...")
				QApplication.processEvents()
				try:
					self._haar_ca, self._haar_cd = haar_dwt(self.seq)
				except Exception as ex:
					QMessageBox.warning(self.main_win, "Wavelet error", str(ex))
					c.btn_temp.setChecked(True); return None
				finally:
					self._sb.showMessage("")
			n_half = self._haar_ca.shape[2]
			img    = self._haar_ca[:, :, min(idx, n_half - 1)].astype(float)

		# ---- Raw temperature ------------------------------------------------
		else:
			img = self.seq[:, :, idx].astype(float)

		img = apply_filter(img, c)
		if crop:
			is_full = (r[0] == 0 and r[2] == 0
					   and r[1] == self.W - 1 and r[3] == self.H - 1)
			if not is_full:
				img = img[r[2]:r[3] + 1, r[0]:r[1] + 1]
		return img

	def _build_lut(self):
		return _build_lut(self.ctrl.cmap_name, self.ctrl.n_colors)

	def _resolve_clamp(self, img_crop, img_full):
		mode = self.ctrl.clamp_mode
		if mode==CLAMP_IMG:
			lo,hi = _nanminmax(img_crop)
		elif mode==CLAMP_SEQ:
			lo,hi = _nanminmax(self.seq) if self.seq is not None else (0.,1.)
		elif mode==CLAMP_ROI:
			r=self.roi; sub=img_full[r[2]:r[3]+1, r[0]:r[1]+1]
			lo,hi = _nanminmax(sub) if sub.size>0 else _nanminmax(img_full)
		else:
			lo,hi=self.ctrl.get_clamp_edits()
			if lo is None: return 0.,1.
			return lo,hi
		if lo==hi: hi=lo+0.1
		self.ctrl.set_clamp(lo,hi)
		return lo,hi

	def _roi_offset(self):
		r=self.roi
		if r[0]==0 and r[2]==0 and r[1]==self.W-1 and r[3]==self.H-1: return None
		return (r[0],r[2])

	def _refresh(self, _=None):
		if self._loading: return
		idx = self._current_idx()
		img_full = self._compute_frame(idx, crop=False)
		if img_full is None: return
		img_crop = self._compute_frame(idx, crop=True)
		if img_crop is None: return
		lut = self._build_lut()
		lo,hi = self._resolve_clamp(img_crop, img_full)
		self.canvas.show_frame(img_crop, (lo, hi), lut, roi_offset=self._roi_offset())
		self.profile_win.update_frame_marker(idx)
		self._push_roi_preview_data(img_full, lo, hi)

	def _push_roi_preview(self, lo=None, hi=None):
		if self.seq is None: return
		idx=int(np.clip(self._current_idx(),0,self.N-1))
		img_full=self._compute_frame(idx,crop=False)
		if img_full is None: return
		if lo is None:
			lo,hi=self.ctrl.get_clamp_edits()
			if lo is None: lo,hi=_nanminmax(img_full)
		self._push_roi_preview_data(img_full,lo,hi)

	def _push_roi_preview_data(self, img_full, lo, hi):
		c=self.ctrl
		d = self.display_panel
		self.roi_canvas.update_view(img_full, self.roi, (lo, hi),
									c.cmap_name,
									d.chk_inv_x.isChecked(),
									d.chk_inv_y.isChecked())

	# -- clamp ----------------------------------------------------------------

	def _on_clamp_mode(self, _=None):
		if self._loading: return
		self.ctrl.set_manual_enabled(self.ctrl.clamp_mode==CLAMP_MAN)
		self._refresh()

	def _do_fit_img(self):
		img=self._compute_frame(self._current_idx(),crop=True)
		if img is None: return
		lo,hi=_nanminmax(img)
		if lo==hi: hi=lo+0.1
		self.ctrl.set_clamp(lo,hi); self.canvas.set_clim(lo,hi)

	def _on_clamp_slider(self):
		if not self.ctrl.btn_clamp_man.isChecked(): return
		lo,hi=self.ctrl.get_clamp_sliders()
		if lo>=hi: lo=hi-1e-9
		self.ctrl.edt_lo.setText(f"{lo:.4g}"); self.ctrl.edt_hi.setText(f"{hi:.4g}")
		self.canvas.set_clim(lo,hi); self._push_roi_preview(lo,hi)

	def _apply_clamp_edits(self):
		lo,hi=self.ctrl.get_clamp_edits()
		if lo is not None:
			if lo==hi: hi=lo+0.1
			self.ctrl.set_clamp(lo,hi); self.canvas.set_clim(lo,hi)
			self._push_roi_preview(lo,hi)

	# -- colour ---------------------------------------------------------------

	def _on_cmap(self):
		idx=self.ctrl.cb_cmap.currentIndex(); cmax=_CMAP_MAXCOLS[idx]
		n=min(self.ctrl.n_colors,cmax)
		self.ctrl.lbl_ncolors.setText(f"Colors: {n} / {cmax}")
		self._refresh()

	def _on_ncolors(self):
		idx=self.ctrl.cb_cmap.currentIndex(); cmax=_CMAP_MAXCOLS[idx]
		n=min(self.ctrl.n_colors,cmax)
		self.ctrl.lbl_ncolors.setText(f"Colors: {n} / {cmax}")
		self._refresh()

	# -- filter ---------------------------------------------------------------

	def _on_filter(self,_=None):
		gauss=self.ctrl.btn_gauss.isChecked()
		self.ctrl.sl_var.setEnabled(gauss); self.ctrl.lbl_var.setEnabled(gauss)
		if not self._loading: self._refresh()

	def _on_variance(self):
		v=self.ctrl.variance; self.ctrl.lbl_var.setText(f"Variance: {v:.2f}")
		if self.ctrl.btn_gauss.isChecked() and not self._loading: self._refresh()

	# -- display --------------------------------------------------------------

	def _on_axis(self, on):
		self.canvas.toggle_axis(on)
		self.display_panel.chk_grid.setEnabled(on)
		if not on: self.display_panel.chk_grid.setChecked(False)

	def _on_frame(self,val):
		self.frame_panel.set_frame_label(val,self.N)
		if not self._loading: self._refresh()

	# -- status bar -----------------------------------------------------------

	def _on_mouse_pos(self,x,y):
		if self.seq is None: return
		off = self._roi_offset()
		ox,oy = off if off else (0,0)
		fc=int(np.clip(round(x+ox),0,self.W-1))
		fr=int(np.clip(round(y+oy),0,self.H-1))
		idx=self._current_idx()
		z=self.seq[fr,fc,idx]
		self._hover_x=fc; self._hover_y=fr; self._hover_z=float(z)
		self._update_idle_status()

	def _on_mouse_out(self):
		self._hover_x = self._hover_y = None
		# _hover_z intentionally kept — shows last known value
		self._update_idle_status()

	# -- canvas click ---------------------------------------------------------

	def _toggle_setupct(self,on):
		self.profile_win.setVisible(on)
		if not on: self.ctrl.btn_extrap.setChecked(False)

	# -- ROI ------------------------------------------------------------------

	def _start_roi(self):
		self._cancel_roi(silent=True)
		self.roi=[0,self.W-1,0,self.H-1]
		self._fft_cache=self._fft_type=self._pct_U=self._pct_S=None
		self._roi_selecting=True; self._roi_pt1=None
		self.canvas.roi_clear_overlay()
		self._sb.showMessage("Click first corner on the main image...")

	def _on_canvas_click(self,x,y):
		if self._roi_selecting:
			if self._roi_pt1 is None:
				self._roi_pt1=(x,y)
				self.canvas.roi_show_pt1(x,y)
				self._sb.showMessage("Click second corner...")
			else:
				self._finish_roi(x,y)
			return
		if not self.ctrl.btn_setupct.isChecked(): return
		off=self._roi_offset(); ox,oy=off if off else (0,0)
		col=int(np.clip(round(x+ox),0,self.W-1))
		row=int(np.clip(round(y+oy),0,self.H-1))
		self.profile_win.add_curve((col,row),self.seq[row,col,:].astype(float))

	def mouseMoveForROI(self,x,y):
		"""Called from canvas_mouse_pos when ROI is active."""
		if not self._roi_selecting or self._roi_pt1 is None: return
		img=self._compute_frame(self._current_idx(),crop=True)
		if img is None: return
		ih,iw=img.shape
		x0,y0=self._roi_pt1
		self.canvas.roi_update_preview(x0,y0,x,y,ih,iw)

	def _finish_roi(self,x1,y1):
		x0,y0=self._roi_pt1; self._roi_selecting=False; self._roi_pt1=None
		self.canvas.roi_clear_overlay()
		off=self._roi_offset(); ox,oy=off if off else (0,0)
		rx0=int(np.clip(round(min(x0,x1)+ox),0,self.W-1))
		rx1=int(np.clip(round(max(x0,x1)+ox),0,self.W-1))
		ry0=int(np.clip(round(min(y0,y1)+oy),0,self.H-1))
		ry1=int(np.clip(round(max(y0,y1)+oy),0,self.H-1))
		if rx1<=rx0 or ry1<=ry0:
			self._sb.showMessage("ROI too small - cancelled."); return
		self.roi=[rx0,rx1,ry0,ry1]
		self._fft_cache=self._fft_type=self._pct_U=self._pct_S=None
		self._sb.showMessage(f"ROI: x [{rx0}:{rx1}]  y [{ry0}:{ry1}]")
		self._refresh()

	def _cancel_roi(self, silent=False):
		self._roi_selecting=False; self._roi_pt1=None
		self.canvas.roi_clear_overlay()
		if not silent: self._sb.showMessage("ROI selection cancelled.")

	def _reset_roi(self):
		self._cancel_roi(silent=True)
		self.roi=[0,self.W-1,0,self.H-1]
		self._fft_cache=self._fft_type=self._pct_U=self._pct_S=None
		self._refresh()

	# -- FFT ------------------------------------------------------------------

	def _get_fft(self,fft_type):
		if self._fft_type==fft_type and self._fft_cache is not None: return self._fft_cache
		self._sb.showMessage("Computing FFT..."); QApplication.processEvents()
		try:
			f=np.fft.fft(self.seq.astype(np.float32),axis=2)
			self._fft_cache=(np.abs(f) if fft_type==1 else -np.angle(f)).astype(np.float32)
			self._fft_type=fft_type
		except MemoryError:
			QMessageBox.warning(self.main_win,"Memory","Not enough memory for FFT."); return None
		finally: self._sb.showMessage("")
		return self._fft_cache

	# -- capture --------------------------------------------------------------

	def _start_capture_dialog(self):
		if self._mode != MODE_STREAM:
			QMessageBox.information(self.main_win,"Capture",
				"Capture is only available in streaming mode."); return
		dlg = CaptureDialog(self.main_win)
		if dlg.exec_() != QDialog.Accepted: return
		path,_=QFileDialog.getSaveFileName(self.main_win,"Save capture","","MATLAB (*.mat)")
		if not path: return
		if not path.endswith(".mat"): path+=".mat"
		n   = dlg.n_samples
		tot = dlg.total_time
		interval_ms = int(tot/n*1000)
		self._cap_buf=[]; self._cap_n=n; self._cap_path=path
		self._cap_timer=QTimer(); self._cap_timer.setInterval(interval_ms)
		self._cap_timer.timeout.connect(self._capture_tick)
		self._cap_timer.start()
		self._sb.showMessage(f"Capturing {n} frames over {tot:.1f}s...")

	def _capture_tick(self):
		if self.seq is not None:
			self._cap_buf.append(self.seq[:,:,0].copy())
		if len(self._cap_buf)>=self._cap_n:
			self._cap_timer.stop(); self._cap_timer=None
			cube=np.stack(self._cap_buf,axis=2).astype(np.float32)
			try:
				save_capture_mat(self._cap_path,cube)
				self._sb.showMessage(
					f"Capture saved: {os.path.basename(self._cap_path)}  "
					f"({cube.shape[0]}x{cube.shape[1]}x{cube.shape[2]})")
			except Exception as e:
				QMessageBox.critical(self.main_win,"Capture error",str(e))

	# -- exports --------------------------------------------------------------

	def _export_png(self):
		if self.seq is None: return
		path,_=QFileDialog.getSaveFileName(self.main_win,"Export PNG","","PNG (*.png)")
		if not path: return
		if not path.endswith(".png"): path+=".png"
		# Grab displayed QPixmap from IRDisplay
		px=self.canvas.grab()
		px.save(path,"PNG")
		self._sb.showMessage(f"Saved {os.path.basename(path)}")

	# -- session --------------------------------------------------------------

	def _collect_state(self):
		c=self.ctrl; pw=self.profile_win
		lo,hi=c.get_clamp_edits()
		mode = next((i for i, b in enumerate(
			[c.btn_temp, c.btn_ftamp, c.btn_ftphase, c.btn_extrap,
			 c.btn_corr, c.btn_ca, c.btn_dac,
			 c.btn_rx, c.btn_skq, c.btn_pca_m, c.btn_tsr, c.btn_wavelet])
					 if b.isChecked()), 0)
		flt=0 if c.btn_noflt.isChecked() else 1 if c.btn_gauss.isChecked() else 2
		return {"filter":flt,"variance":c.variance,"cmap":c.cb_cmap.currentIndex(),
				"ncolors":c.n_colors,"clamp_lo":lo or 0.,"clamp_hi":hi or 1.,
				"clamp_mode":c.clamp_mode,"inv_x":c.chk_inv_x.isChecked(),
				"inv_y":c.chk_inv_y.isChecked(),"axis":c.chk_axis.isChecked(),
				"grid":c.chk_grid.isChecked(),"proc":mode,
				"frame":self._current_idx(),"roi":self.roi,
				"et":pw.etype,"ig":pw.sl_ig.value(),"if_":pw.sl_if.value(),
				"dg":pw.sl_dg.value(),"df":pw.sl_df.value(),
				"cg":pw.sl_cg.value(),"cf":pw.sl_cf.value(),
				"Ls":pw.sl_L.value(),"hs":pw.sl_h.value()}

	def _save_session(self):
		if self.seq is None:
			QMessageBox.warning(self.main_win,"Save","No data loaded."); return
		path,_=QFileDialog.getSaveFileName(self.main_win,"Save session","","IR View session (*.irv)")
		if not path: return
		if not path.endswith(".irv"): path+=".irv"
		self._sb.showMessage("Saving...")
		try: save_irv(path,self._collect_state(),self.seq)
		except Exception as e: QMessageBox.critical(self.main_win,"Save error",str(e))
		finally: self._sb.showMessage("")
		self.main_win.setWindowTitle(f"IR View - {os.path.basename(path)}")

	def _load_session(self,path):
		try: state,seq=load_irv(path)
		except Exception as e:
			QMessageBox.critical(self.main_win,"Load error",str(e)); return
		self._source_path=path; self._mode=MODE_FILE
		self._ingest(seq)
		self._loading=True; c=self.ctrl; pw=self.profile_win
		try:
			self.roi=list(state.get("roi",[0,self.W-1,0,self.H-1]))
			[c.btn_noflt,c.btn_gauss,c.btn_median][state.get("filter",0)].setChecked(True)
			c.sl_var.setValue(int(state.get("variance",0.85)*100))
			c.cb_cmap.setCurrentIndex(state.get("cmap",9))
			c.sl_ncolors.setValue(state.get("ncolors",256))
			c.set_clamp(state.get("clamp_lo",0.),state.get("clamp_hi",1.))
			c.chk_inv_x.setChecked(state.get("inv_x",False))
			c.chk_inv_y.setChecked(state.get("inv_y",True))
			c.chk_axis.setChecked(state.get("axis",True))
			c.chk_grid.setChecked(state.get("grid",False))
			[c.btn_temp, c.btn_ftamp, c.btn_ftphase, c.btn_extrap,
			 c.btn_corr, c.btn_ca, c.btn_dac,
			 c.btn_rx, c.btn_skq, c.btn_pca_m, c.btn_tsr,
			 c.btn_wavelet][state.get("proc", 0)].setChecked(True)
			fp=self.frame_panel
			fp.sl_frame.setValue(int(np.clip(state.get("frame",0),0,self.N-1))+1)
			cmode=state.get("clamp_mode",CLAMP_IMG)
			[c.btn_clamp_img,c.btn_clamp_seq,c.btn_clamp_roi,c.btn_clamp_man][cmode].setChecked(True)
			c.set_manual_enabled(cmode==CLAMP_MAN)
			pw.sl_ig.setValue(state.get("ig",0));  pw.sl_if.setValue(state.get("if_",0))
			pw.sl_dg.setValue(state.get("dg",0));  pw.sl_df.setValue(state.get("df",10))
			pw.sl_cg.setValue(state.get("cg",0));  pw.sl_cf.setValue(state.get("cf",20))
			pw.sl_L.setValue(state.get("Ls",400)); pw.sl_h.setValue(state.get("hs",100))
			pw.cb_etype.setCurrentIndex(state.get("et",1)-1)
		finally: self._loading=False
		self.main_win.setWindowTitle(f"IR View - {os.path.basename(path)}")
		self._refresh()

	# -- entry point ----------------------------------------------------------

	def run(self):
		# Wire mouse-move to ROI rubber-band
		self.canvas.canvas_mouse_pos.connect(self.mouseMoveForROI)
		self.main_win.show()
		sys.exit(self.app.exec_())

# ---------------------------------------------------------------------------

def main():
	iv=IrView()
	if len(sys.argv)>1: iv._load_mat(sys.argv[1])
	iv.run()

if __name__=="__main__":
	main()
