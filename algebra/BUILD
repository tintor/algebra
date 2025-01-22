package(default_visibility = ["//visibility:public"])
load(":library.bzl", "library")

cc_library(name = "__test", hdrs = ["__test.h"], deps = ["@catch2//:catch2_main"])

library("natural",       hdrs = ["integer_backend.h"])
library("natural_func",  deps = [":natural"])
library("integer",       deps = [":natural"])
library("integer_func",  deps = [":integer"])
library("rational",      deps = [":integer_func"])
library("rational_func", deps = [":rational", ":natural_func"])
library("real",          deps = [":rational", ":natural_func"])
library("real_func",     deps = [":real"])
library("expr",          deps = [":rational_func"])

library("vector",          has_test = False)
library("rational_vector", has_test = False, deps = [":vector", ":rational_func"])
library("solve_linear",    has_test = False, deps = [":vector"])

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
    deps = [":natural_func", ":integer_func", ":rational_func", ":real_func", ":expr", ":geometry", ":rational_vector"],
)
