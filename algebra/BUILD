package(default_visibility = ["//visibility:public"])

cc_library(name = "integer_backend", hdrs = ["integer_backend.h"])
cc_library(name = "natural", hdrs = ["natural.h"], deps = [":integer_backend"])
cc_library(name = "integer", hdrs = ["integer.h"], deps = [":natural"])
cc_library(name = "rational", hdrs = ["rational.h"], deps = [":integer"])
cc_library(name = "real", hdrs = ["real.h"], deps = [":rational"])
cc_library(name = "expr", hdrs = ["expr.h"], deps = [":rational"])

cc_test(
    name = "natural_test",
    srcs = ["natural_test.cc"],
    deps = [":natural", "@catch2//:catch2_main"],
    args = ["-d=yes"],
)

cc_test(
    name = "integer_test",
    srcs = ["integer_test.cc"],
    deps = [":integer", "@catch2//:catch2_main"],
    args = ["-d=yes"],
)

cc_test(
    name = "rational_test",
    srcs = ["rational_test.cc"],
    deps = [":rational", "@catch2//:catch2_main"],
    args = ["-d=yes"],
)

cc_test(
    name = "real_test",
    srcs = ["real_test.cc"],
    deps = [":real", "@catch2//:catch2_main"],
    args = ["-d=yes"],
)

cc_test(
    name = "expr_test",
    srcs = ["expr_test.cc"],
    deps = [":expr", "@catch2//:catch2_main"],
    args = ["-d=yes"],
)
