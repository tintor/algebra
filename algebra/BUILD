package(default_visibility = ["//visibility:public"])
load(":library.bzl", "library")

cc_library(name = "__test", hdrs = ["__test.h"], deps = ["@catch2//:catch2_main", ":kernels"])

library("kernels",       hdrs = ["util.h", "types.h"])
library("natural_class", hdrs = ["integer_backend.h"], deps = [":kernels"])
library("natural",       deps = [":natural_class"])
library("integer_class", deps = [":natural"])
library("integer",       deps = [":integer_class"])
library("rational_class",deps = [":integer"])
library("rational",      deps = [":rational_class"])
library("xrational",     deps = [":rational"])
library("real_class",    deps = [":rational"])
library("real",          deps = [":real_class"])
library("expr",          deps = [":rational"])

library("vector",          test_deps = [":rational"])
library("solve_linear",    deps = [":vector"], test_deps = [":rational"])
library("rational_vector", deps = [":vector", ":rational", ":xrational"])

library("point_segment_squared_distance",
    deps = [":vector"],
    test_deps = [":rational_vector"])

library("segment_segment_squared_distance",
    deps = [":solve_linear"],
    test_deps = [":rational_vector"])

library("segment_segment_intersection",
    deps = [":solve_linear"],
    test_deps = [":rational_vector"])

library("polygon2",
    deps = [":vector", ":rational"],
    test_deps = [":rational_vector"])

library("polygon2_arc",
    deps = [":polygon2", ":solve_linear"],
    test_deps = [":rational_vector"])

library("polygon2_boolean",
    deps = [":polygon2", ":segment_segment_intersection"],
    test_deps = [":rational_vector"])

library("polygon2_buffer",
    deps = [":polygon2_boolean"],
    test_deps = [":rational_vector"])

library("polygon2_arc_boolean",
    deps = [":polygon2_arc"],
    test_deps = [":rational_vector"])

library("geometry",
    deps = [":solve_linear", ":point_segment_squared_distance", ":segment_segment_squared_distance", ":segment_segment_intersection"],
    test_deps = [":rational_vector"])

library("dual")

cc_library(
    name = "algebra",
    hdrs = ["algebra.h"],
    deps = [":rational", ":real", ":expr", ":geometry", ":polygon2", ":polygon2_arc", ":polygon2_arc_boolean", ":polygon2_boolean", ":polygon2_buffer", ":rational_vector", ":xrational", ":dual"],
)

cc_test(
    name = "integer_backend_test",
    srcs = ["integer_backend_test.cc"],
    deps = [":natural_class", ":__test"],
    args = ["-d=yes"],
)

cc_test(
    name = "portable_divq_test",
    srcs = ["portable_divq_test.cc"],
    deps = [":kernels", ":__test"],
    args = ["-d=yes"],
)

cc_test(
    name = "link_test",
    srcs = ["link_test.cc", "link_test2.cc"],
    deps = [":algebra", ":__test"],
    args = ["-d=yes"],
)

cc_binary(name = "long_stress_test", srcs = ["long_stress_test.cc", "__stress_test.h"], deps = [":rational", ":xrational", ":real"])
cc_binary(name = "isqrt_benchmark", srcs = ["isqrt_benchmark.cc"], deps = ["@google_benchmark//:benchmark", ":natural"])
