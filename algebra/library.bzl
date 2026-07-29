# One header, its library, and the tests that go with it.
#
#   library("integer") declares cc_library integer from integer.h, and cc_test integer_test from
#   integer_test.cc against it.
#
# tests: extra test files against the same library, for a header whose tests are split over several
#   files (integer_class has three). Each entry `x` is a cc_test named x from x.cc.
# srcs: extra source files for a test, for one that needs a second translation unit (link_test).
def library(name, hdrs=[], srcs=[], deps=[], test_deps=[], has_test=True, tests=[], test_srcs={}):
    native.cc_library(
        name = name,
        hdrs = [name + ".h"] + hdrs,
        srcs = srcs,
        deps = deps,
    )

    all_tests = ([name + "_test"] if has_test else []) + tests
    for test in all_tests:
        native.cc_test(
            name = test,
            srcs = [test + ".cc"] + test_srcs.get(test, []),
            deps = test_deps + [":" + name, ":__test"],
            args = ["-d=yes"],
        )
