#include "stencil.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <functional>
#include <immintrin.h>
#include <iostream>
#include <mutex>
#include <random>
#include <thread>
#include <utility>
#include <vector>

std::mutex mu[20000];
std::thread threads[10];

// 400ms

template <typename P>
void calc(P *in, P *out, int height, int width, int from, int to) {
    for (int i = from; i <= to; i++) {
        for (int j = 1; j < width - 1; j++) {
            int im1jm1 = (i - 1) * width + j - 1;
            int im1j = (i - 1) * width + j;
            int im1jp1 = (i - 1) * width + j + 1;
            int ijm1 = (i)*width + j - 1;
            int ij = (i)*width + j;
            int ijp1 = (i)*width + j + 1;
            int ip1jm1 = (i + 1) * width + j - 1;
            int ip1j = (i + 1) * width + j;
            int ip1jp1 = (i + 1) * width + j + 1;
            P val = -in[im1jm1] - in[im1j] - in[im1jp1] - in[ijm1] + 8 * in[ij] - in[ijp1] - in[ip1jm1] - in[ip1j] -
                    in[ip1jp1];
            val = (val < 0 ? 0 : val);
            val = (val > 255 ? 255 : val);
            out[i * width + j] = val;
        }
    }
}

float tmp[] = {8.0f, 8.0f, 8.0f, 8.0f};

// 3500ms...
/*
template <typename P>
void calc(P *in, P *out, int height, int width, int from, int to) {
    for (int i = 0; i < 9; i++) {
        int po = from * width + 1, pi;
        switch (i) {
            case 0:
                pi = (from - 1) * width;
                break;
            case 1:
                pi = (from - 1) * width + 1;
                break;
            case 2:
                pi = (from - 1) * width + 2;
                break;
            case 3:
                pi = (from)*width;
                break;
            case 4:
                pi = (from)*width + 2;
                break;
            case 5:
                pi = (from + 1) * width;
                break;
            case 6:
                pi = (from + 1) * width + 1;
                break;
            case 7:
                pi = (from + 1) * width + 2;
                break;
            case 8:  // sepcial case
                pi = (from)*width;
        }

        for (int k = 0; k <= to - from; k++) {
            int p, q;
            for (p = pi + k * width, q = po + k * width; p < pi + k * width + width - 2; p += 4, q += 4) {
                if (i < 8) {
                    __m128 u = _mm_loadu_ps(in + p);
                    __m128 v = _mm_loadu_ps(out + q);
                    __m128 r = _mm_sub_ps(v, u);
                    _mm_storeu_ps(out + q, r);
                } else {
                    __m128 u = _mm_loadu_ps(in + p);
                    __m128 t = _mm_loadu_ps(tmp);
                    u = _mm_mul_ps(u, t);
                    __m128 v = _mm_loadu_ps(out + q);
                    __m128 r = _mm_add_ps(v, u);
                    _mm_storeu_ps(out + q, r);
                }
            }
        }
    }
}
*/

/*
// might be overflow
template <typename P>
void calc(P *in, P *out, std::mutex mu[], const int height, const int width, int id) {
    int po = (1 + id * height / 9) * width + 1, pi;
    switch (id) {
        case 0:
            pi = (1 + id * height / 9 - 1) * width;
            break;
        case 1:
            pi = (1 + id * height / 9 - 1) * width + 1;
            break;
        case 2:
            pi = (1 + id * height / 9 - 1) * width + 2;
            break;
        case 3:
            pi = (1 + id * height / 9) * width;
            break;
        case 4:
            pi = (1 + id * height / 9) * width + 2;
            break;
        case 5:
            pi = (1 + id * height / 9 + 1) * width;
            break;
        case 6:
            pi = (1 + id * height / 9 + 1) * width + 1;
            break;
        case 7:
            pi = (1 + id * height / 9 + 1) * width + 2;
            break;
        case 8:  // sepcial case
            pi = (1 + id * height / 9) * width + 1;
    }

    // with or without lock, all ~300ms in my machine?
    for (int i = 1; i < height - 1; i++) {
        // mu[i].lock();

        for (int j = 0; j < width - 2; j++) {
            if (id < 8)
                out[po] -= in[pi];
            else
                out[po] += in[pi] * 8;
            po++;
            pi++;
        }

        po += 2;
        pi += 3;
        if (po >= height * width)
            po -= height * width;
        if (pi >= height * width)
            pi -= height * width;
        // mu[i].unlock();
    }
}
*/

template <typename P>
void ApplyStencil(ImageClass<P> &img_in, ImageClass<P> &img_out) {
    const int width = img_in.width;
    const int height = img_in.height;

    P *in = img_in.pixel;
    P *out = img_out.pixel;

    int thread_num = 8;
    /*
    thread_num = 9;
    for (int i = 0; i < thread_num; i++) {
        threads[i] = std::thread(calc<P>, in, out, mu, height, width, i);
    }
    */

    // for (int i = 0; i < 9; i++)
    // calc<P>(in, out, mu, height, width, i);

    for (int i = 0; i < thread_num; i++)
        threads[i] = std::thread(calc<P>,
                                 in,
                                 out,
                                 height,
                                 width,
                                 (height - 2) / thread_num * i + 1,
                                 i == thread_num - 1 ? height - 2 : (height - 2) / thread_num * (i + 1));

    for (int i = 0; i < thread_num; i++)
        threads[i].join();
}

template void ApplyStencil<float>(ImageClass<float> &img_in, ImageClass<float> &img_out);