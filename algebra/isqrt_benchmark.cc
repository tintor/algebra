#include "algebra/natural.h"
#include <benchmark/benchmark.h>
using namespace benchmark;

#include <random>
#include <print>

using std::format;
using std::print;
using namespace algebra;
using namespace algebra::literals;

using ucent = unsigned __int128;
using cent = __int128;
using ulong = unsigned long;
using uint = unsigned int;

namespace std {
inline namespace __1 {
bool __is_posix_terminal(__sFILE*) { return true; }
}
}

static void test(const auto& fn, benchmark::State& state) {
    std::mt19937_64 rng(0);
    const int bits = state.range(0);
    std::vector<natural> x;
    const int N = 16;
    for (int i = 0; i < N; i++)
        x.push_back(uniform_sample_bits(bits, rng));

    int i = 0;
    for (auto _ : state) {
        DoNotOptimize(fn(x[i]));
        i = ++i % N;
    }
}

#define BENCH_LOW(FUNCTION) \
static void BM_LOW_ ## FUNCTION(benchmark::State& state) { test(FUNCTION, state); } \
BENCHMARK(BM_LOW_ ## FUNCTION)->RangeMultiplier(2)->Range(16, 131072);

BENCH_LOW(isqrt2)
BENCH_LOW(isqrt_simple)
BENCH_LOW(isqrt_newthon)
BENCH_LOW(isqrt_binary_search)
BENCH_LOW(isqrt_digits)

BENCHMARK_MAIN();
