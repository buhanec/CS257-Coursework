https://www.cs.fsu.edu/~engelen/courses/HPC-adv/MMXandSSEexamples.txt

================================================================================

Problem: vectorize the following code with MMX and SSE

char a[N], b[N], c[N];
...
for (i = 0; i < N; i++)
  a[i] = b[i] + c[i];

MMX technology with 64-bit MM registers (aligned load/store):

Back:
  movq  mm0, _b[ecx]
  paddb mm0, _c[ecx]
  movq  _a[ecx], mm0
  add   ecx, 8
  cmp   ecx, edi
  jl    Back

MMX intrinsics (aligned load/store):

char a[N], b[N], c[N];
...
__m64 *av, *bv, *cv;
av = (__m64*)a; // assume 8-byte aligned
bv = (__m64*)b; // assume 8-byte aligned
cv = (__m64*)c; // assume 8-byte aligned
for (i = 0; i < N/8; i++)
  av[i] = _mm_add_pi8(bv[i], cv[i]);

SSE/SSE2 technology with 128-bit XMM registers (aligned load/store):

Back:
  movdqa xmm0, _b[ecx]
  paddb  xmm0, _c[ecx]
  movdqa _a[ecx], xmm0
  add    ecx, 16
  cmp    ecx, edi
  jl     Back

SSE2 intrinsics (aligned load/store):

char a[N], b[N], c[N];
...
__m128i *av, *bv, *cv;
av = (__m128i*)a; // assume 16-byte aligned
bv = (__m128i*)b; // assume 16-byte aligned
cv = (__m128i*)c; // assume 16-byte aligned
for (i = 0; i < N/16; i++)
  av[i] = _mm_add_epi8(bv[i], cv[i]);

SSE/SSE2 technology with 128-bit XMM registers (unaligned load/store):

Back:
  movdqu xmm0, _b[ecx]
  movdqu xmm1, _c[ecx]
  paddb  xmm0, xmm1
  movdqu _a[ecx], xmm0
  add    ecx, 16
  cmp    ecx, edi
  jl     Back

SSE2 intrinsics (unaligned load/store)

char a[N], b[N], c[N];
...
__m128i *av, *bv, *cv;
av = (__m128i*)a;
bv = (__m128i*)b;
cv = (__m128i*)c;
for (i = 0; i < N/16; i++)
{ __m128i br = _mm_loadu_si128(&bv[i];
  __m128i cr = _mm_loadu_si128(&cv[i];
  __m128i ar = _mm_add_epi8(br, cr);
  _mm_storeu_si128(&av[i], ar);
}

================================================================================

Problem: vectorize the following code with SSE

double a[N], b[N], c[N];
...
for (i = 0; i < N; i++)
  a[i] = b[i] + c[i];

SSE/SSE2 technology with 128-bit XMM registers (aligned load/store):

Back:
  movapd xmm0, _b[ecx]      |b1   |b0   |
  addpd  xmm0, _c[ecx]      |b1+c1|b0+c0|
  movapd _a[ecx], xmm0
  add    ecx, 16
  cmp    ecx, edi
  jl     Back

SSE2 intrinsics (aligned load/store):

double a[N], b[N], c[N];
...
__m128d *av, *bv, *cv;
av = (__m128d*)a; // assume 16-byte aligned
bv = (__m128d*)b; // assume 16-byte aligned
cv = (__m128d*)c; // assume 16-byte aligned
for (i = 0; i < N/2; i++)
  av[i] = _mm_add_pd(bv[i], cv[i]);

================================================================================

Problem: vectorize the following code with SSE

float a[N], b[N], x;
...
for (i = 0; i < N; i++)
  a[i] = b[i] + x;

SSE/SSE2 technology with 128-bit XMM registers (aligned load/store):

  movss  xmm0, _x       |0|0|0|x|
  shufps xmm0, xmm0, 0      |x|x|x|x|
Back:
  movaps xmm1, _b[ecx]      |b3  |b2  |b1  |b0  |
  addps  xmm1, xmm0     |b3+x|b2+x|b1+x|b0+x|
  movaps _a[ecx], xmm1
  add    ecx, 16
  cmp    ecx, edi
  jl     Back

Note: instead of shufps we can twice unpack:

  movss    xmm0, _x     |0|0|0|x|
  unpcklps xmm0, xmm0       |0|0|x|x|
  unpcklps xmm0, xmm0       |x|x|x|x|

SSE2 intrinsics (aligned load/store):

float a[N], b[N], x;
...
__m128 *av, *bv, xr;
av = (__m128d*)a; // assume 16-byte aligned
bv = (__m128d*)b; // assume 16-byte aligned
xr = _mm_set1_ps(x);
for (i = 0; i < N/4; i++)
  av[i] = _mm_add_ps(bv[i], xr);

================================================================================

Problem: scalar expand the following code with SSE

double x;
int n;

SSE/SSE2 technology with 128-bit XMM registers (aligned load/store):

  movsd    xmm0, _x     |0|x|
  unpcklpd xmm0, xmm0       |x|x|
  movdqa   xmm1, _n     |0|0|0|n|
  unpcklwd xmm1, xmm1       |0|0|n|n|
  unpcklwd xmm1, xmm1       |n|n|n|n|

SSE2 intrinsics (aligned load/store):

__m128d xr = _mm_set1_pd(x);
__m128i nr = _mm_set1_epi32(n);

================================================================================

Problem: vectorize the following code with SSE

char a[N];
...
for (i = 0; i < N; i++)
  a[i] = i;

SSE/SSE2 technology with 128-bit XMM registers (aligned load/store):

  movdqa xmm0, _cnst$1      |15|14|13|12|11|10| 9| 8| 7| 6| 5| 4| 3| 2| 1| 0|
  movdqa xmm1, _cnst$2      |16|16|16|16|16|16|16|16|16|16|16|16|16|16|16|16|
Back:
  movdqa _a[ecx], xmm0
  paddb  xmm0, xmm1
  add    ecx, 16
  cmp    ecx, edi
  jl     Back

SSE2 intrinsics (aligned load/store):

char a[N];
...
__m128i *av = (__m128i*)a; // assume 16-byte aligned
__m128i iv = _mm_set_epi8(15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0)
__m128i ci = _mm_set1_epi8(16);
av = (__m128i*)a; // assume 16-byte aligned
for (i = 0; i < N/16; i++)
{ av[i] = iv;
  iv = _mm_add_epi8(iv, ci);
}

================================================================================

Problem: vectorize the following code with SSE

double a[N], x = 0.0;
...
for (i = 0; i < N; i++)
{
  x = x + 4;
  a[i] = 6 * x + 1;
}

SSE/SSE2 technology with 128-bit XMM registers (aligned load/store):

  movapd xmm0, _cnst$1      |49.0|25.0|
  movapd xmm1, _cnst$2      |48.0|48.0|
Back:
  movapd _a[ecx], xmm0
  addpd  xmm0, xmm1
  add    ecx, 16
  cmp    ecx, edi
  jl     Back

SSE2 intrinsics (aligned load/store):

double a[N], x = 0.0;
...
__m128d *av = (__m128d*)a; // assume 16-byte aligned
__m128d iv = _mm_set_pd(49.0,25.0)
__m128d ci = _mm_set1_pd(48.0);
av = (__m128d*)a; // assume 16-byte aligned
for (i = 0; i < N/2; i++)
  av[i] = iv = _mm_add_epi8(iv, ci);

================================================================================

Problem: vectorize the following code with SSE

int a[N], x = 0;
...
for (i = 0; i < N; i++)
  x = x + a[i];

SSE/SSE2 technology with 128-bit XMM registers (aligned load/store):

  pxor    xmm0, xmm0        |0|0|0|0|
Back:
  movdqa  xmm1, _a[ecx]
  paddd   xmm0, xmm1
  add     ecx, 16
  cmp     ecx, edi
  jl      Back
  psrldq  xmm1, xmm0, 32    |0  |0    |x3      |x2         |
  paddw   xmm0, xmm1            |x3 |x2   |x1+x3   |x0+x2      |
  psrldq  xmm1, xmm0, 16        |0  |x3   |x2      |x1+x3      |
  paddw   xmm0, xmm1            |x3 |x2+x3|x1+x2+x3|x0+x1+x2+x3|


int a[N], x = 0;
...
SSE2 intrinsics (aligned load/store):

int a[N], x = 0.0;
...
int xx[4];
__m128i *av = (__m128i*)a; // assume 16-byte aligned
__m128i xv = _mm_setzero_si128();
for (i = 0; i < N/4; i++)
  xv = _mm_add_epi32(xv, av[i]);
_mm_store_si128(xx, xv);
x = xx[0] + xx[1] + xx[2] + xx[3];

================================================================================

Problem: vectorize the following code with SSE

float a[N], b[N], x = 0.0;
...
for (i = 0; i < N; i++)
  x = x + a[i]*b[i];

SSE/SSE2 technology with 128-bit XMM registers (aligned load/store):

  xorps  xmm0, xmm0     |0|0|0|0|
Back:
  movaps xmm1, _a[ecx]
  mulps  xmm1, _b[ecx]
  addps  xmm0, xmm1
  add    ecx, 16
  cmp    ecx, edi
  jl     Back
  haddps xmm0, xmm0     |x3+x2      |x1+x0      |x3+x2      |x2+x0      |
  haddps xmm0, xmm0     |x3+x2+x1+x0|x3+x2+x1+x0|x3+x2+x1+x0|x3+x2+x1+x0|

SSE2 intrinsics (aligned load/store):

float a[N], x = 0.0;
...
float xx[4];
__m128 *av = (__m128*)a; // assume 16-byte aligned
__m128 *bv = (__m128*)b; // assume 16-byte aligned
__m128 xv = _mm_setzero_ps();
for (i = 0; i < N/4; i++)
  xv = _mm_add_ps(xv, _mm_mul_ps(av[i], bv[i]));
_mm_store_ps(xx, xv);
x = xx[0] + xx[1] + xx[2] + xx[3];
// or instead of the two lines above we can use a horizontal add:
xv = _mm_hadd_ps(xv, xv);
xv = _mm_hadd_ps(xv, xv);
_mm_store_ps(xx, xv);
x = xx[0];

================================================================================

Problem: vectorize the following code with SSE

double a[N], x = 0.0;
...
for (i = 0; i < N; i++)
  x = x + a[3*i];

SSE/SSE2 technology with 128-bit XMM registers (unaligned load/store):

  xorpd  xmm0, xmm0     |0|0|
Back:
  movsd  xmm1, _a[ecx]      |0      |a[3*i+0]|
  movhpd xmm1, _a[ecx+24]   |a[3*i+3|a[3*i+0]|
  addpd  xmm0, xmm1
  add    ecx, 48
  cmp    ecx, edi
  jl     Back
  haddpd xmm0, xmm0     |x1+x0|x1+x0|

SSE2 intrinsics (unaligned load/store):

double a[N], x = 0.0;
...
double xx[2];
__m128d xv = _mm_setzero_pd();
__m128d t;
for (i = 0; i < N; i += 2)
{
  t = _mm_load_sd(a+3*i);         // need not be aligned
  t = _mm_loadh_pd(t, a+3*(i+1)); // need not be aligned
  xv = xv + t;
}
_mm_store_pd(xx, xv);
x = xx[0] + xx[1];

================================================================================

Problem: vectorize the following code with SSE

float a[N], x = 0.0;
...
for (i = 0; i < N; i++)
  x = x + a[3*i];

SSE/SSE2 technology with 128-bit XMM registers (unaligned load/store):

  xorps    xmm0, xmm0       |0|0|0|0|
Back:
  movss    xmm4, _a[ecx]    |-       |-       |-       |a[3*i+0]|
  movss    xmm3, _a[ecx+12] |-       |-       |-       |a[3*i+3]|
  movss    xmm2, _a[ecx+24] |-       |-       |-       |a[3*i+6]|
  movss    xmm1, _a[ecx+36] |-       |-       |-       |a[3*i+9]|
  unpcklps xmm4, xmm2       |-       |-       |a[3*i+6]|a[3*i+0]|
  unpcklps xmm3, xmm1       |-       |-       |a[3*i+9]|a[3*i+3]|
  unpcklps xmm4, xmm3       |a[3*i+9]|a[3*i+6]|a[3*i+3]|a[3*i+0]|
  addps    xmm0, xmm4
  add      ecx, 48
  cmp      ecx, edi
  jl       Back
  haddps   xmm0, xmm0       |x3+x2      |x1+x0      |x3+x2      |x2+x0      |
  haddps   xmm0, xmm0       |x3+x2+x1+x0|x3+x2+x1+x0|x3+x2+x1+x0|x3+x2+x1+x0|

SSE2 intrinsics (unaligned load/store):

float a[N], x = 0.0;
...
float xx[4];
__m128 xv = _mm_setzero_ps();
__m128 t;
for (i = 0; i < N; i += 4)
{
  t = _mm_set_ps(a+3*i, a+3*i+3, a+3*i+6, a+3*i+9);
  xv = _mm_add_ps(xv, t);
}
_mm_store_ps(xx, xv);
x = xx[0] + xx[1] + xx[2] + xx[3];

Note: performance may be poor, since this is not much better than scalar code!

================================================================================

Problem: vectorize the following code with SSE

float a[N], b[N], c[N];
...
for (i = 0; i < N; i++)
  if (a[i] > 0)
    a[i] = b[i] / c[i];

SSE/SSE2 technology with 128-bit XMM registers (aligned load/store):

Back:
  movaps   xmm0, _a[ecx]
  movaps   xmm2, _b[ecx]
  divps    xmm2, _c[ecx]    set  x[] = |b3/c3   |b2/c2   |b1/c1   |b0/c0   |
  xorps    xmm1, xmm1       set 0
  cmpltps  xmm1, xmm0       guards g[]=|a3>0    |a2>0    |a1>0    |a0>0    |
  movaps   xmm3, xmm1       copy guards g[]
  andnps   xmm3, xmm0       mask y[] = |!a3>0&a3|!a2>0&a2|!a1>0&a1|!a0>0&a0|
  andps    xmm2, xmm1       mask z[] = | a3>0&x3| a2>0&x2| a1>0&x1| a0>0&x0|
  orps     xmm3, xmm2       combine  = |y3|z3   |y2|z2   |y1|z1   |y0|z0   |
  movaps   _a[ecx], xmm3    store into a[]
  add      ecx, 16
  cmp      ecx, edi
  jl       Back

SSE2 intrinsics (aligned load/store):

float a[N], b[N], c[N];
...
__m128 *av = (__m128*)a; // assume 16-byte aligned
__m128 *bv = (__m128*)b; // assume 16-byte aligned
__m128 *cv = (__m128*)c; // assume 16-byte aligned
__m128 zeros = _mm_setzero_ps();
for (i = 0; i < N/4; i++)
{
  __m128 x = _mm_div_ps(bv[i], cv[i]);
  __m128 g = _mm_cmplt_ps(av[i], zeros);
  __m128 y = _mm_andnot_ps(g, av[i]);
  __m128 z = _mm_and_ps(g, x);
  av[i] = _mm_or_ps(y, z);
}
