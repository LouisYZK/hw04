#include <cstdio>
#include <cstdlib>
#include <vector>
#include <chrono>
#include <cmath>

float frand() {
    return (float)rand() / RAND_MAX * 2 - 1;
}

struct StarVec {
    std::vector<float> px, py, pz;
    std::vector<float> vx, vy, vz;
    std::vector<float> mass;
};

StarVec stars;

void init() {
    for (int i = 0; i < 48; i++) {
        stars.px.push_back(frand());
        stars.py.push_back(frand());
        stars.pz.push_back(frand());
        stars.vx.push_back(frand());
        stars.vy.push_back(frand());
        stars.vz.push_back(frand());
        stars.mass.push_back(frand() + 1);
    }
}

float G = 0.001;
float eps = 0.001;
float dt = 0.01;

void step() {
    for ( size_t i = 0; i < 48; ++i ) {
        float cum_dx = 0, cum_dy = 0, cum_dz = 0;
        #pragma opm simd
        for ( size_t j = 0; j < 48; ++j ) {
            float dx = stars.px[j] - stars.px[i];
            float dy = stars.py[j] - stars.py[i];
            float dz = stars.pz[j] - stars.pz[i];
            float d2 = dx * dx + dy * dy + dz * dz + eps * eps;
            d2 *= std::sqrt(d2);
            cum_dx+= dx * stars.mass[j] * G * dt / d2;
            cum_dy += dy * stars.mass[j] * G * dt / d2;
            cum_dz += dz * stars.mass[j] * G * dt / d2;
        }
        stars.vx[i] += cum_dx;
        stars.vy[i] += cum_dy;
        stars.vz[i] += cum_dz;
    }
    for ( size_t i = 0; i < 48; ++i ) {
        stars.px[i] += stars.vx[i] * dt;
        stars.py[i] += stars.vy[i] * dt;
        stars.pz[i] += stars.vz[i] * dt;
    }
}

float calc() {
    float energy = 0;
    for ( size_t i = 0; i < 48; ++i ) {
        float v2 = stars.vx[i] * stars.vx[i] + stars.vy[i] * stars.vy[i] + stars.vz[i] * stars.vz[i];
        energy += stars.mass[i] * v2 / 2;
        for ( size_t j = 0; j < 48; ++j ) {
            float dx = stars.px[j] - stars.px[i];
            float dy = stars.py[j] - stars.py[i];
            float dz = stars.pz[j] - stars.pz[i];
            float d2 = dx * dx + dy * dy + dz * dz + eps * eps;
            energy -= stars.mass[j] * stars.mass[i] * G / sqrt(d2) / 2;
        }
    }
    return energy;
}

template <class Func>
long benchmark(Func const &func) {
    auto t0 = std::chrono::steady_clock::now();
    func();
    auto t1 = std::chrono::steady_clock::now();
    auto dt = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0);
    return dt.count();
}

int main() {
    init();
    printf("Initial energy: %f\n", calc());
    auto dt = benchmark([&] {
        for (int i = 0; i < 100000; i++)
            step();
    });
    printf("Final energy: %f\n", calc());
    printf("Time elapsed: %ld ms\n", dt);
    return 0;
}
