package(default_visibility = ["//visibility:public"])
load(":library.bzl", "library")

cc_library(name = "__test", hdrs = ["__test.h"], deps = ["@catch2//:catch2_main"])

library("kernels",       hdrs = ["util.h", "types.h"])
library("natural_class", hdrs = ["integer_backend.h"], deps = [":kernels"])
library("natural",       deps = [":natural_class"])
library("integer_class", deps = [":natural"])
library("integer",       deps = [":integer_class"])
library("rational",      deps = [":integer"])
library("xrational",     deps = [":rational"])
library("rational_func", deps = [":rational"])
library("real_class",    deps = [":rational"])
library("real",          deps = [":real_class"])
library("expr",          deps = [":rational_func"])

library("vector",          has_test = False)
library("solve_linear",    has_test = False, deps = [":vector"])
library("rational_vector", deps = [":vector", ":rational_func", ":xrational"])

library("point_segment_squared_distance",
    deps = [":vector"],
    test_deps = [":rational_vector"])

library("segment_segment_squared_distance",
    deps = [":solve_linear"],
    test_deps = [":rational_vector"])

library("segment_segment_intersection",
    deps = [":solve_linear"],
    test_deps = [":rational_vector"])

library("geometry",
    deps = [":solve_linear", ":point_segment_squared_distance", ":segment_segment_squared_distance", ":segment_segment_intersection"],
    test_deps = [":rational_vector"])

cc_library(
    name = "algebra",
    deps = [":integer", ":rational_func", ":real", ":expr", ":geometry", ":rational_vector"],
)

cc_binary(name = "long_stress_test", srcs = ["long_stress_test.cc", "__stress_test.h"], deps = [":integer", ":rational_func", ":xrational", ":real"])
cc_binary(name = "isqrt_benchmark", srcs = ["isqrt_benchmark.cc"], deps = ["@google_benchmark//:benchmark", ":natural"])
