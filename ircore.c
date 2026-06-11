/*
 *  This file is part of IR View Suite.
 *
 *  IR View Suite - Thermal/IR imaging suite
 *  Copyright (C) 2026 Hugo Barriga
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License v3
 *  as published by the Free Software Foundation.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY.
 *  See the GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU AGPLv3 along with this program.
 *  If not, see <https://www.gnu.org/licenses/>.
 */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <float.h>

static inline float _at(const float *a, int r, int c, int H, int W)
{
	int rr = (r < 0) ? 0 : (r >= H ? H - 1 : r);
	int cc = (c < 0) ? 0 : (c >= W ? W - 1 : c);
	return a[rr * W + cc];
}

static void _isort9(float *a)
{
	#define SW(i, j) do { if (a[i] > a[j]) { float t = a[i]; a[i] = a[j]; a[j] = t; } } while (0)

	SW(0, 1); SW(3, 4); SW(6, 7);
	SW(1, 2); SW(4, 5); SW(7, 8);

	SW(0, 1); SW(3, 4); SW(6, 7);
	SW(0, 3); SW(3, 6); SW(0, 3);

	SW(1, 4); SW(4, 7); SW(1, 4);
	SW(2, 5); SW(5, 8); SW(2, 5);

	SW(1, 3); SW(5, 7); SW(2, 6);
	SW(4, 6); SW(2, 4); SW(2, 3);
	SW(4, 5); SW(5, 6); SW(3, 4);

	#undef SW
}

static inline float _median9(float *v)
{
	_isort9(v);
	return v[4];
}

static inline float _median3(float a, float b, float c)
{
	if (a > b) { float t = a; a = b; b = t; }
	if (b > c) { float t = b; b = c; c = t; }
	if (a > b) { float t = a; a = b; b = t; }
	return b;
}

void hmf(const float *src, float *dst, int H, int W)
{
	static const int CR[9][2] = {
		{ -2,  0 }, { -1, 0 }, { 0, -2 },
		{  0, -1 }, {  0, 0 }, { 0,  1 },
		{  0,  2 }, {  1, 0 }, { 2,  0 }
	};

	static const int DI[9][2] = {
		{ -2, -2 }, { -1, -1 }, {  0,  0 },
		{  1,  1 }, {  2,  2 }, { -2,  2 },
		{ -1,  1 }, {  1, -1 }, {  2, -2 }
	};

	float *tmp = (float *)malloc((size_t)H * W * sizeof(float));
	if (!tmp) {
		memcpy(dst, src, (size_t)H * W * sizeof(float));
		return;
	}

	for (int r = 0; r < H; r++) {
		for (int c = 0; c < W; c++) {

			float cr[9], di[9];

			for (int k = 0; k < 9; k++) {
				cr[k] = _at(src, r + CR[k][0], c + CR[k][1], H, W);
				di[k] = _at(src, r + DI[k][0], c + DI[k][1], H, W);
			}

			tmp[r * W + c] = _median3(_median9(cr), _median9(di), src[r * W + c]);
		}
	}

	memcpy(dst, tmp, (size_t)H * W * sizeof(float));
	free(tmp);
}

void gauss(const float *src, float *dst, int H, int W, float variance)
{
	double i_max = 10.0 * sqrt((double)variance);
	int klen = (int)round(i_max) + 1;

	float *kern = (float *)malloc((size_t)klen * sizeof(float));
	if (!kern) {
		memcpy(dst, src, (size_t)H * W * sizeof(float));
		return;
	}

	double half = i_max * 0.5;
	double inv2v = 1.0 / (2.0 * (double)variance);
	double norm = 1.0 / sqrt(2.0 * M_PI * (double)variance);

	double sum = 0.0;
	for (int i = 0; i < klen; i++) {
		double d = i - half;
		double v = norm * exp(-d * d * inv2v);
		kern[i] = (float)v;
		sum += v;
	}

	for (int i = 0; i < klen; i++) {
		kern[i] = (float)(kern[i] / sum);
	}

	int bw = klen / 2;
	int off = klen - 1 - bw;

	int PH = H + 2 * bw;
	int PW = W + 2 * bw;

	float *pad = (float *)calloc((size_t)PH * PW, sizeof(float));
	float *tmp2 = (float *)calloc((size_t)PH * PW, sizeof(float));

	if (!pad || !tmp2) {
		free(kern);
		free(pad);
		free(tmp2);
		memcpy(dst, src, (size_t)H * W * sizeof(float));
		return;
	}

	// Replicate padding
	for (int r = 0; r < H; r++) {
		int pr = r + bw;

		for (int c = 0; c < bw; c++) {
			pad[pr * PW + c] = src[r * W];
		}

		memcpy(&pad[pr * PW + bw], &src[r * W], (size_t)W * sizeof(float));

		for (int c = 0; c < bw; c++) {
			pad[pr * PW + bw + W + c] = src[r * W + W - 1];
		}
	}

	for (int pr = 0; pr < bw; pr++) {
		memcpy(&pad[pr * PW], &pad[bw * PW], (size_t)PW * sizeof(float));
	}

	for (int pr = bw + H; pr < PH; pr++) {
		memcpy(&pad[pr * PW], &pad[(bw + H - 1) * PW], (size_t)PW * sizeof(float));
	}

	// Horizontal pass
	for (int r = 0; r < PH; r++) {
		for (int c = 0; c < PW; c++) {

			double acc = 0.0;

			const float *row = &pad[r * PW];

			for (int k = 0; k < klen; k++) {
				int sc = c + off - k;
				if (sc < 0) sc = 0;
				else if (sc >= PW) sc = PW - 1;

				acc += row[sc] * kern[k];
			}

			tmp2[r * PW + c] = (float)acc;
		}
	}

	// Vertical pass
	for (int r = 0; r < H; r++) {
		int rr = r + bw;

		for (int c = 0; c < W; c++) {

			double acc = 0.0;

			for (int k = 0; k < klen; k++) {
				int sr = rr + off - k;
				if (sr < 0) sr = 0;
				else if (sr >= PH) sr = PH - 1;

				acc += tmp2[sr * PW + (c + bw)] * kern[k];
			}

			dst[r * W + c] = (float)acc;
		}
	}

	free(kern);
	free(pad);
	free(tmp2);
}

void corrct(const float *seq, float *out, int H, int W, int N)
{
	int HW = H * W;

	double *ref = (double *)calloc((size_t)N, sizeof(double));
	if (!ref) {
		memset(out, 0, (size_t)HW * sizeof(float));
		return;
	}

	for (int t = 0; t < N; t++) {
		double s = 0.0;
		const float *col = &seq[t];
		for (int p = 0; p < HW; p++) {
			s += col[p * N];
		}
		ref[t] = s / (double)HW;
	}

	double ref_mean = 0.0;
	for (int t = 0; t < N; t++) ref_mean += ref[t];
	ref_mean /= (double)N;

	double ref_var = 0.0;
	for (int t = 0; t < N; t++) {
		double d = ref[t] - ref_mean;
		ref_var += d * d;
	}

	double ref_norm = sqrt(ref_var);

	for (int p = 0; p < HW; p++) {
		const float *row = &seq[p * N];

		double pm = 0.0;
		for (int t = 0; t < N; t++) pm += row[t];
		pm /= (double)N;

		double num = 0.0, pv = 0.0;

		for (int t = 0; t < N; t++) {
			double pc = row[t] - pm;
			double rc = ref[t] - ref_mean;
			num += pc * rc;
			pv += pc * pc;
		}

		double den = sqrt(pv) * ref_norm;
		double corr = (den > 1e-300) ? (num / den) : 0.0;

		double val = 1.0 - corr;
		double v5 = val * val * val * val * val;

		out[p] = (float)cbrt(v5);
	}

	free(ref);
}

void pct_standardize(float *mat, int rows, int cols)
{
	for (int c = 0; c < cols; c++) {

		double mu = 0.0;

		for (int r = 0; r < rows; r++) {
			mu += mat[r * cols + c];
		}

		mu /= (double)rows;

		double var = 0.0;

		for (int r = 0; r < rows; r++) {
			double d = mat[r * cols + c] - mu;
			var += d * d;
		}

		double sig = (rows > 1) ? sqrt(var / (double)(rows - 1)) : 1.0;
		if (sig < 1e-300) sig = 1.0;

		float fmu = (float)mu;
		float fsig = (float)sig;

		for (int r = 0; r < rows; r++) {
			float *v = &mat[r * cols + c];
			*v = (*v - fmu) / fsig;
		}
	}
}

void extrap_sib(const float *img_prime, double t_p, double t_now, float *out, int n)
{
	float ratio = (t_now > 0.0) ? (float)sqrt(t_p / t_now) : 1.0f;

	for (int i = 0; i < n; i++) {
		out[i] = img_prime[i] * ratio;
	}
}

void minmax(const float *a, int n, float *out_min, float *out_max)
{
	float lo = FLT_MAX;
	float hi = -FLT_MAX;

	for (int i = 0; i < n; i++) {
		float v = a[i];
		if (v != v) continue;

		if (v < lo) lo = v;
		if (v > hi) hi = v;
	}

	*out_min = lo;
	*out_max = hi;
}

void ca(const float *seq, const float *cold, int ref_px, float *out, int HW, int N)
{
	for (int i = 0; i < N; i++) {

		float ref_val = seq[ref_px * N + i] - cold[ref_px];

		for (int p = 0; p < HW; p++) {
			out[p * N + i] = (seq[p * N + i] - cold[p]) - ref_val;
		}
	}
}

void dac(const float *seq, const float *cold, int t_prime, float *out, int HW, int N)
{
	int out_N = N - t_prime;
	if (out_N <= 0) return;

	int tp0 = t_prime - 1;

	for (int i = 0; i < out_N; i++) {

		double ratio = sqrt((double)t_prime / (double)(i + 1));

		for (int p = 0; p < HW; p++) {

			float Bpi  = seq[p * N + (i + tp0)] - cold[p];
			float Bref = seq[p * N + tp0] - cold[p];

			out[p * out_N + i] = (float)(Bpi - ratio * Bref);
		}
	}
}

void rx(const float *seq_centred, const float *sphere, float *out, int HW, int N)
{
	for (int p = 0; p < HW; p++) {

		double score = 0.0;

		const float *x = &seq_centred[p * N];

		for (int k = 0; k < N; k++) {

			double yk = 0.0;

			const float *sk = &sphere[k * N];

			for (int j = 0; j < N; j++) {
				yk += sk[j] * x[j];
			}

			score += yk * yk;
		}

		out[p] = (float)score;
	}
}

void skq(const float *seq_centred, float *kurt, float *skew, float *fifth, int HW, int N)
{
	for (int p = 0; p < HW; p++) {

		const float *x = &seq_centred[p * N];

		double m2 = 0.0, m3 = 0.0, m4 = 0.0, m5 = 0.0;

		for (int i = 0; i < N; i++) {
			double v = x[i];
			double v2 = v * v;

			m2 += v2;
			m3 += v2 * v;
			m4 += v2 * v2;
			m5 += v2 * v2 * v;
		}

		m2 /= N; m3 /= N; m4 /= N; m5 /= N;

		double sig2 = m2;
		double sig = sqrt(sig2);

		double sig3 = sig2 * sig;
		double sig4 = sig2 * sig2;

		kurt[p] = (sig4 > 1e-30) ? (float)(m4 / sig4) : 0.0f;
		skew[p] = (sig3 > 1e-30) ? (float)(m3 / sig3) : 0.0f;

		double sig5 = sig4 * sig;
		fifth[p] = (sig5 > 1e-30) ? (float)(m5 / sig5) : 0.0f;
	}
}

void haar_dwt(const float *seq, float *ca, float *cd, int HW, int N)
{
	int half = N / 2;
	const float S = 0.7071067811865476f;

	for (int p = 0; p < HW; p++) {

		const float *x = &seq[p * N];
		float *a = &ca[p * half];
		float *d = &cd[p * half];

		for (int k = 0; k < half; k++) {
			float x0 = x[2 * k];
			float x1 = x[2 * k + 1];

			a[k] = S * (x0 + x1);
			d[k] = S * (x0 - x1);
		}
	}
}
