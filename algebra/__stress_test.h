#include <print>
#include <random>
#include <vector>
#include <string>
#include <thread>
#include <atomic>
#include <mutex>
#include <source_location>
#include <any>
#include <cxxabi.h>
using std::print;
using std::format;

void reset_color() {
    print("\033[1;0m");
}

void set_red() {
    print("\033[1;31m");
}

void set_blue() {
    print("\033[1;34m");
}

void set_yellow() {
    print("\033[1;33m");
}

std::mutex g_mutex;
int g_remaining_errors = 0;

#define TEST(E) { M __m{#E, seed, std::source_location::current(), []() -> std::string {return "";}}; catch_exceptions(__m, [&](){__m <=> E;}); }
#define TEST2(E, MSG) { M __m{#E, seed, std::source_location::current(), [&]() -> std::string { return MSG;}}; catch_exceptions(__m, [&](){__m <=> E;}); }

template<typename Msg>
struct M {
    const char* expr;
    uint64_t seed;
    std::source_location loc;
    Msg msg;

    void print_failure() const {
        print("{}:{} with seed {} ", loc.file_name(), loc.line(), seed);
        set_red();
        print("FAILED:\n");
        set_blue();
        print("  TEST( {} )\n", expr);
        std::string m = msg();
        if (m.size()) {
            reset_color();
            print("with message:\n");
            set_yellow();
            std::stringstream ss(m);
            std::string line;
            while (getline(ss, line, '\n'))
                print("  {}\n", line);
        }
        reset_color();
    }
};

struct handled_error { };

void try_report_error(const auto& fn) {
    {
        std::lock_guard g(g_mutex);
        fn();
        if (--g_remaining_errors <= 0)
            exit(0);
    }
    throw handled_error();
}

template<typename T>
void catch_exceptions(const M<T>& m, const auto& fn) {
    try {
        fn();
    } catch (handled_error& e) {
        throw;
    } catch (std::runtime_error& e) {
        try_report_error([&](){
            m.print_failure();
            print("due to unexpected exception with message:\n");
            print("  {}\n", e.what());
        });
    } catch (...) {
        try_report_error([&](){
            m.print_failure();
            print("due to unexpected exception\n");
        });
    }
}

std::string shorten(const std::string& a, int pre, int post) {
    return (a.size() <= pre + post) ? a : (a.substr(0, pre) + "..." + a.substr(a.size() - post));
}

template<typename A, typename T>
struct M1 {
    M<T> m;
    A a;

    template<typename B>
    void print_failure(const char* op, const B& b) const {
        m.print_failure();
        print("with expansion:\n");
        set_yellow();
        auto as = shorten(std::format("{}", a), 100, 100);
        auto bs = shorten(std::format("{}", b), 100, 100);
        if (as.size() >= 20 || bs.size() >= 20)
            print("  {}\n  {}\n  {}\n", as, op, bs);
        else
            print("  {} {} {}\n", as, op, bs);
        reset_color();
    }
};

template<typename A, typename T>
M1<A, T> operator<=>(const M<T>& m, const A& a) { return {m, a}; }

template<typename T>
void operator<=>(const M<T>& m, const bool value) {
    if (!value) {
        try_report_error([&](){
            std::lock_guard g(g_mutex);
            m.print_failure();
        });
    }
}

template<typename A, typename T>
void catch_exceptions(const M1<A, T>& m1, const char* op, const auto& fn, const auto& b) {
    try {
        if (!fn())
            try_report_error([&](){
                m1.print_failure(op, b);
            });
    } catch (handled_error& e) {
        throw;
    } catch (std::runtime_error& e) {
        try_report_error([&](){
            m1.print_failure(op, b);
            print("due to unexpected exception with message:\n");
            print("  {}\n", e.what());
        });
    } catch (...) {
        try_report_error([&](){
            m1.print_failure(op, b);
            print("due to unexpected exception\n");
        });
    }
}

#define M_OP(OP) template<typename A, typename T> void operator OP (const M1<A, T>& m1, const auto& b) { catch_exceptions(m1, #OP, [&](){ return m1.a OP b; }, b); }

M_OP(==)
M_OP(!=)
M_OP(<)
M_OP(>)
M_OP(<=)
M_OP(>=)

void stress_test(const auto& fn, int max_errors = 6) {
    std::random_device rd;
    std::atomic<uint64_t> seed = (uint64_t(rd()) << 32) + rd();
    g_remaining_errors = max_errors;

    auto func = [&](){
        while (true) {
            uint64_t s = seed++;
            try {
                try {
                    fn(seed);
                } catch (handled_error&) {
                } catch (std::runtime_error& e) {
                    try_report_error([&](){
                        reset_color();
                        print("UKNOWN with seed {} ", s);
                        set_red();
                        print("FAILED:\n");
                        reset_color();
                        print("with message:\n");
                        set_yellow();
                        print("  {}\n", e.what());
                        reset_color();
                    });
                } catch (...) {
                    int status;
                    char* demangled = abi::__cxa_demangle(abi::__cxa_current_exception_type()->name(), 0, 0, &status);

                    try_report_error([&](){
                        reset_color();
                        print("UKNOWN with seed {} ", s);
                        set_red();
                        print("FAILED:\n");
                        reset_color();
                        if (demangled) {
                            print("with exception_type:\n");
                            set_yellow();
                            print("  {}\n", demangled);
                            reset_color();
                        }
                    });
                    free(demangled);
                }
            } catch (handled_error&) { }
            if (s % 1000 == 0) {
                std::lock_guard g(g_mutex);
                print("seed {}\n", s);
            }
        }
    };

    std::vector<std::thread> threads;
    for (int i = 0; i < std::thread::hardware_concurrency(); i++)
        threads.push_back(std::thread(func));
    threads[0].join();
}

namespace std {
inline namespace __1 {
bool __is_posix_terminal(__sFILE*) { return true; }
}
}
