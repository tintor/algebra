package(default_visibility = ["//visibility:public"])

cc_library(name = "integer_backend", hdrs = ["integer_backend.h"])
cc_library(name = "natural", hdrs = ["natural.h"], deps = [":integer_backend"])
cc_library(name = "integer", hdrs = ["integer.h"], deps = [":natural"])
cc_library(name = "rational", hdrs = ["rational.h"], deps = [":integer"])
cc_library(name = "real", hdrs = ["real.h"], deps = [":rational"])
cc_library(name = "expr", hdrs = ["expr.h"], deps = [":rational"])

cc_library(name = "natural_func", hdrs = ["natural_func.h"], deps = [":natural"])
cc_library(name = "integer_func", hdrs = ["integer_func.h"], deps = [":integer"])
cc_library(name = "rational_func", hdrs = ["rational_func.h"], deps = [":rational"])
cc_library(name = "real_func", hdrs = ["real_func.h"], deps = [":real"])

cc_library(name = "__test", hdrs = ["__test.h"], deps = ["@catch2//:catch2_main"])

cc_test(
    name = "natural_test",
    srcs = ["natural_test.cc"],
    deps = [":natural", ":__test"],
    args = ["-d=yes"],
)

cc_test(
    name = "integer_test",
    srcs = ["integer_test.cc"],
    deps = [":integer", ":__test"],
    args = ["-d=yes"],
)

cc_test(
    name = "rational_test",
    srcs = ["rational_test.cc"],
    deps = [":rational", ":__test"],
    args = ["-d=yes"],
)

cc_test(
    name = "real_test",
    srcs = ["real_test.cc"],
    deps = [":real", ":__test"],
    args = ["-d=yes"],
)

cc_test(
    name = "expr_test",
    srcs = ["expr_test.cc"],
    deps = [":expr", ":__test"],
    args = ["-d=yes"],
)

cc_test(
    name = "natural_func_test",
    srcs = ["natural_func_test.cc"],
    deps = [":natural_func", ":__test"],
    args = ["-d=yes"],
)

cc_test(
    name = "integer_func_test",
    srcs = ["integer_func_test.cc"],
    deps = [":integer_func", ":__test"],
    args = ["-d=yes"],
)

cc_test(
    name = "rational_func_test",
    srcs = ["rational_func_test.cc"],
    deps = [":rational_func", ":__test"],
    args = ["-d=yes"],
)

cc_test(
    name = "real_func_test",
    srcs = ["real_func_test.cc"],
    deps = [":real_func", ":__test"],
    args = ["-d=yes"],
)
