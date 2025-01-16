#pragma once
#include <catch2/catch_test_macros.hpp>
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
