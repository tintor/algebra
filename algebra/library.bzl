def library(name, hdrs=[], srcs=[], deps=[], test_deps=[], has_test=True):
    native.cc_library(
        name = name,
        hdrs = [name + ".h"] + hdrs,
        srcs = srcs,
        deps = deps,
    )
    if has_test:
        native.cc_test(
            name = name + "_test",
            srcs = [name + "_test.cc"],
            deps = test_deps + [":" + name, ":__test"],
            args = ["-d=yes"],
        )
