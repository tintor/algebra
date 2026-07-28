#include "algebra/integer_backend.h"
#include "algebra/__test.h"

// integer_backend is a small buffer optimization: it holds either a single inline word
// (allocated() == false, capacity() == 1) or a heap buffer, with the sign folded into the
// size. Every test below says which side of that boundary it is on, since that is where
// these operations can go wrong.

TEST_CASE("integer_backend - default state") {
    integer_backend a;
    REQUIRE(a.size() == 0);
    REQUIRE(a.sign() == 0);
    REQUIRE(a.empty());
    REQUIRE(!a.allocated());
    REQUIRE(a.capacity() == 1);
    REQUIRE(a[0] == 0);
    REQUIRE(a.data()[0] == 0);
}

TEST_CASE("integer_backend - a single word stays inline") {
    integer_backend a {5};
    REQUIRE(!a.allocated());
    REQUIRE(a.capacity() == 1);
    REQUIRE(a.size() == 1);
    REQUIRE(a.sign() == 1);
    REQUIRE(a[0] == 5);

    integer_backend b(-5);
    REQUIRE(!b.allocated());
    REQUIRE(b.size() == 1);
    REQUIRE(b.sign() == -1);
    REQUIRE(b[0] == 5);
}

TEST_CASE("integer_backend - data points at the words") {
    integer_backend a {4, 5, 6};
    REQUIRE(a.allocated());
    REQUIRE(a.data()[0] == 4);
    REQUIRE(a.data()[1] == 5);
    REQUIRE(a.data()[2] == 6);
    a.data()[1] = 9;
    REQUIRE(a[1] == 9);
    REQUIRE(static_cast<const integer_backend&>(a).data()[1] == 9);

    integer_backend b {7};
    REQUIRE(!b.allocated());
    REQUIRE(b.data()[0] == 7);
    b.data()[0] = 8;
    REQUIRE(b[0] == 8);
    REQUIRE(static_cast<const integer_backend&>(b).data()[0] == 8);
}

TEST_CASE("integer_backend - capacity of an inline value is one") {
    integer_backend a;
    REQUIRE(a.capacity() == 1);
    REQUIRE(!a.allocated());
    a.push_back(7);
    REQUIRE(a.capacity() == 1);
    REQUIRE(!a.allocated());
    a.push_back(8); // no longer fits inline
    REQUIRE(a.allocated());
    REQUIRE(a.capacity() >= 2);
    REQUIRE(a.size() == 2);
    REQUIRE(a[0] == 7);
    REQUIRE(a[1] == 8);
}

TEST_CASE("integer_backend - capacity is never smaller than the size") {
    integer_backend a;
    for (uint64_t i = 0; i < 40; i++) {
        a.push_back(i);
        REQUIRE(a.capacity() >= a.size());
        REQUIRE(a.size() == static_cast<int>(i) + 1);
    }
    for (int i = 0; i < 40; i++)
        REQUIRE(a[i] == static_cast<uint64_t>(i));
}

// ----- resize -----

TEST_CASE("integer_backend - resize to one word stays inline") {
    integer_backend a;
    a.resize(1);
    REQUIRE(a.size() == 1);
    REQUIRE(a.sign() == 1);
    REQUIRE(!a.allocated());
    REQUIRE(a.capacity() == 1);
    REQUIRE(a[0] == 0);
}

TEST_CASE("integer_backend - resize clears a stale inline word") {
    integer_backend a {5};
    a.pop_back(); // size 0, but the inline word still holds 5
    REQUIRE(a.size() == 0);
    a.resize(1);
    REQUIRE(a.size() == 1);
    REQUIRE(a[0] == 0);
}

TEST_CASE("integer_backend - resize from inline to heap") {
    integer_backend a {7};
    a.resize(3);
    REQUIRE(a.allocated());
    REQUIRE(a.capacity() >= 3);
    REQUIRE(a.size() == 3);
    REQUIRE(a.sign() == 3);
    REQUIRE(a[0] == 7); // the inline word becomes word 0
    REQUIRE(a[1] == 0);
    REQUIRE(a[2] == 0);
}

TEST_CASE("integer_backend - resize keeps the sign") {
    integer_backend a(-7);
    a.resize(3);
    REQUIRE(a.sign() == -3);
    REQUIRE(a.size() == 3);
    REQUIRE(a[0] == 7);
    a.resize(1);
    REQUIRE(a.sign() == -1);
    REQUIRE(a[0] == 7);
}

TEST_CASE("integer_backend - resize down and up again zeroes the new words") {
    integer_backend a {1, 2, 3};
    a.resize(1);
    REQUIRE(a.size() == 1);
    REQUIRE(a[0] == 1);
    REQUIRE(a.allocated()); // shrinking does not release the buffer
    a.resize(3);
    REQUIRE(a.size() == 3);
    REQUIRE(a[0] == 1);
    REQUIRE(a[1] == 0);
    REQUIRE(a[2] == 0);
}

TEST_CASE("integer_backend - resize grows beyond the capacity") {
    integer_backend a {1, 2};
    a.resize(9);
    REQUIRE(a.size() == 9);
    REQUIRE(a.capacity() >= 9);
    REQUIRE(a[0] == 1);
    REQUIRE(a[1] == 2);
    for (int i = 2; i < 9; i++)
        REQUIRE(a[i] == 0);
}

TEST_CASE("integer_backend - resize to zero") {
    integer_backend a {1, 2};
    a.resize(0);
    REQUIRE(a.size() == 0);
    REQUIRE(a.empty());
    REQUIRE(a.sign() == 0);

    integer_backend b(-5);
    b.resize(0);
    REQUIRE(b.empty());
    REQUIRE(b.sign() == 0);
}

// ----- downsize -----

TEST_CASE("integer_backend - downsize keeps the words and the buffer") {
    integer_backend a {1, 2, 3};
    const int cap = a.capacity();
    a.downsize(2);
    REQUIRE(a.size() == 2);
    REQUIRE(a.sign() == 2);
    REQUIRE(a.capacity() == cap);
    REQUIRE(a.allocated());
    REQUIRE(a[0] == 1);
    REQUIRE(a[1] == 2);
    REQUIRE(a.back() == 2);
    a.downsize(0);
    REQUIRE(a.empty());
    REQUIRE(a.sign() == 0);
}

TEST_CASE("integer_backend - downsize keeps the sign") {
    integer_backend a {1, 2, 3};
    a.negate();
    a.downsize(2);
    REQUIRE(a.sign() == -2);
    REQUIRE(a.size() == 2);
    a.downsize(1);
    REQUIRE(a.sign() == -1);
    REQUIRE(a[0] == 1);
}

// ----- reserve -----

TEST_CASE("integer_backend - reserve of one word stays inline") {
    integer_backend a {5};
    a.reserve(1);
    REQUIRE(!a.allocated());
    REQUIRE(a.capacity() == 1);
    REQUIRE(a.size() == 1);
    REQUIRE(a.sign() == 1);
    REQUIRE(a[0] == 5);
}

TEST_CASE("integer_backend - reserve moves the inline word to the heap") {
    integer_backend a {5};
    a.reserve(4);
    REQUIRE(a.allocated());
    REQUIRE(a.capacity() >= 4);
    REQUIRE(a.size() == 1); // reserve does not change the size
    REQUIRE(a.sign() == 1);
    REQUIRE(a[0] == 5);
}

TEST_CASE("integer_backend - reserve keeps size, sign and words") {
    integer_backend a {1, 2, 3};
    a.negate();
    a.reserve(10);
    REQUIRE(a.allocated());
    REQUIRE(a.capacity() >= 10);
    REQUIRE(a.size() == 3);
    REQUIRE(a.sign() == -3);
    REQUIRE(a[0] == 1);
    REQUIRE(a[1] == 2);
    REQUIRE(a[2] == 3);
}

TEST_CASE("integer_backend - reserve of a negative inline word") {
    integer_backend a(-5);
    a.reserve(3);
    REQUIRE(a.allocated());
    REQUIRE(a.size() == 1);
    REQUIRE(a.sign() == -1);
    REQUIRE(a[0] == 5);
}

TEST_CASE("integer_backend - reserve of an empty value") {
    integer_backend a;
    a.reserve(4);
    REQUIRE(a.allocated());
    REQUIRE(a.capacity() >= 4);
    REQUIRE(a.size() == 0);
    REQUIRE(a.sign() == 0);
    REQUIRE(a.empty());
}

TEST_CASE("integer_backend - reserve does not shrink") {
    integer_backend a {1, 2, 3, 4};
    const int cap = a.capacity();
    const uint64_t* p = a.data();
    a.reserve(0);
    a.reserve(1);
    a.reserve(2);
    a.reserve(cap);
    REQUIRE(a.capacity() == cap);
    REQUIRE(a.data() == p); // no reallocation
    REQUIRE(a.size() == 4);
    REQUIRE(a[3] == 4);
}

TEST_CASE("integer_backend - reserve_bits") {
    integer_backend a {5};
    a.reserve_bits(0);
    REQUIRE(!a.allocated());
    REQUIRE(a.capacity() == 1);
    a.reserve_bits(64); // one word is enough
    REQUIRE(!a.allocated());
    REQUIRE(a.capacity() == 1);
    a.reserve_bits(65);
    REQUIRE(a.allocated());
    REQUIRE(a.capacity() >= 2);
    REQUIRE(a.size() == 1);
    REQUIRE(a[0] == 5);
    a.reserve_bits(64 * 5);
    REQUIRE(a.capacity() >= 5);
    a.reserve_bits(64 * 5 + 1);
    REQUIRE(a.capacity() >= 6);
    REQUIRE(a.size() == 1);
    REQUIRE(a.sign() == 1);
    REQUIRE(a[0] == 5);
}

// ----- reserve_and_set_zero -----

TEST_CASE("integer_backend - reserve_and_set_zero of one word stays inline") {
    integer_backend a(-5);
    a.reserve_and_set_zero(1);
    REQUIRE(a.empty());
    REQUIRE(a.sign() == 0);
    REQUIRE(a[0] == 0);
    REQUIRE(!a.allocated());
    REQUIRE(a.capacity() == 1);
}

TEST_CASE("integer_backend - reserve_and_set_zero within the existing buffer") {
    integer_backend a {1, 2, 3, 4};
    a.negate();
    const int cap = a.capacity();
    const uint64_t* p = a.data();
    a.reserve_and_set_zero(2);
    REQUIRE(a.empty());
    REQUIRE(a.sign() == 0);
    REQUIRE(a[0] == 0);
    REQUIRE(a.data() == p); // the buffer is reused
    REQUIRE(a.capacity() == cap);
}

TEST_CASE("integer_backend - reserve_and_set_zero grows the buffer") {
    integer_backend a {1, 2};
    a.reserve_and_set_zero(9);
    REQUIRE(a.empty());
    REQUIRE(a.size() == 0);
    REQUIRE(a.sign() == 0);
    REQUIRE(a.allocated());
    REQUIRE(a.capacity() >= 9);
}

// ----- pop_back -----

TEST_CASE("integer_backend - pop_back") {
    integer_backend a {1, 2, 3};
    a.pop_back();
    REQUIRE(a.size() == 2);
    REQUIRE(a.sign() == 2);
    REQUIRE(a.back() == 2);
    a.pop_back();
    REQUIRE(a.size() == 1);
    REQUIRE(a.back() == 1);
    a.pop_back();
    REQUIRE(a.empty());
    REQUIRE(a.sign() == 0);
    a.pop_back(); // does not underflow
    REQUIRE(a.empty());
    REQUIRE(a.sign() == 0);
    REQUIRE(a.size() == 0);
}

TEST_CASE("integer_backend - pop_back keeps the sign") {
    integer_backend a {1, 2, 3};
    a.negate();
    a.pop_back();
    REQUIRE(a.sign() == -2);
    REQUIRE(a.size() == 2);
    a.pop_back();
    REQUIRE(a.sign() == -1);
    a.pop_back();
    REQUIRE(a.sign() == 0); // zero has no sign
    a.pop_back();
    REQUIRE(a.sign() == 0);
}

TEST_CASE("integer_backend - pop_back of the inline word") {
    integer_backend a(-5);
    REQUIRE(!a.allocated());
    a.pop_back();
    REQUIRE(a.empty());
    REQUIRE(a.sign() == 0);
    REQUIRE(!a.allocated());
}

// ----- insert_first_word -----

TEST_CASE("integer_backend - insert_first_word into an empty value") {
    integer_backend a;
    a.insert_first_word(9);
    REQUIRE(a.size() == 1);
    REQUIRE(a.sign() == 1);
    REQUIRE(a[0] == 9);
    REQUIRE(!a.allocated()); // one word still fits inline
}

TEST_CASE("integer_backend - insert_first_word moves the inline word to the heap") {
    integer_backend a {7};
    a.insert_first_word(9);
    REQUIRE(a.allocated());
    REQUIRE(a.size() == 2);
    REQUIRE(a[0] == 9);
    REQUIRE(a[1] == 7);
}

TEST_CASE("integer_backend - insert_first_word grows a full heap buffer") {
    integer_backend a {1, 2, 3};
    while (a.size() < a.capacity())
        a.push_back(4);
    const int n = a.size();
    a.insert_first_word(9);
    REQUIRE(a.size() == n + 1);
    REQUIRE(a.capacity() >= n + 1);
    REQUIRE(a[0] == 9);
    REQUIRE(a[1] == 1);
    REQUIRE(a[2] == 2);
    REQUIRE(a[3] == 3);
}

TEST_CASE("integer_backend - insert_first_word keeps the sign") {
    integer_backend a(-7);
    a.insert_first_word(9);
    REQUIRE(a.sign() == -2);
    REQUIRE(a[0] == 9);
    REQUIRE(a[1] == 7);
    a.insert_first_word(0);
    REQUIRE(a.sign() == -3);
    REQUIRE(a[0] == 0);
    REQUIRE(a[1] == 9);
    REQUIRE(a[2] == 7);
}

// ----- insert_first_n_words -----

TEST_CASE("integer_backend - insert_first_n_words of zero is a no-op") {
    integer_backend a {1, 2};
    a.negate();
    a.insert_first_n_words(0);
    REQUIRE(a.size() == 2);
    REQUIRE(a.sign() == -2);
    REQUIRE(a[0] == 1);
    REQUIRE(a[1] == 2);
}

TEST_CASE("integer_backend - insert_first_n_words shifts the words up") {
    integer_backend a {1, 2};
    a.insert_first_n_words(2);
    REQUIRE(a.size() == 4);
    REQUIRE(a.sign() == 4);
    REQUIRE(a[0] == 0);
    REQUIRE(a[1] == 0);
    REQUIRE(a[2] == 1);
    REQUIRE(a[3] == 2);
}

TEST_CASE("integer_backend - insert_first_n_words from inline to heap") {
    integer_backend a {7};
    a.insert_first_n_words(1);
    REQUIRE(a.allocated());
    REQUIRE(a.size() == 2);
    REQUIRE(a[0] == 0);
    REQUIRE(a[1] == 7);
}

TEST_CASE("integer_backend - insert_first_n_words moves the inline word far up") {
    integer_backend a(-7);
    a.insert_first_n_words(3);
    REQUIRE(a.allocated());
    REQUIRE(a.size() == 4);
    REQUIRE(a.sign() == -4);
    REQUIRE(a[0] == 0);
    REQUIRE(a[1] == 0);
    REQUIRE(a[2] == 0);
    REQUIRE(a[3] == 7);
}

TEST_CASE("integer_backend - insert_first_n_words into an empty value") {
    integer_backend a;
    a.insert_first_n_words(1);
    REQUIRE(a.size() == 1);
    REQUIRE(a[0] == 0);
    REQUIRE(!a.allocated());

    integer_backend b;
    b.insert_first_n_words(3);
    REQUIRE(b.size() == 3);
    REQUIRE(b.allocated());
    REQUIRE(b[0] == 0);
    REQUIRE(b[1] == 0);
    REQUIRE(b[2] == 0);
}

TEST_CASE("integer_backend - insert_first_n_words keeps the sign") {
    integer_backend a {1, 2};
    a.negate();
    a.insert_first_n_words(2);
    REQUIRE(a.sign() == -4);
    REQUIRE(a[0] == 0);
    REQUIRE(a[1] == 0);
    REQUIRE(a[2] == 1);
    REQUIRE(a[3] == 2);
}

// ----- erase_first_n_words -----

TEST_CASE("integer_backend - erase_first_n_words") {
    integer_backend a {1, 2, 3};
    a.erase_first_n_words(0); // no-op
    REQUIRE(a.size() == 3);
    REQUIRE(a[0] == 1);
    a.erase_first_n_words(1);
    REQUIRE(a.size() == 2);
    REQUIRE(a.sign() == 2);
    REQUIRE(a[0] == 2);
    REQUIRE(a[1] == 3);
    a.erase_first_n_words(2);
    REQUIRE(a.empty());
    REQUIRE(a.sign() == 0);
}

TEST_CASE("integer_backend - erase_first_n_words keeps the sign") {
    integer_backend a {1, 2, 3};
    a.negate();
    a.erase_first_n_words(2);
    REQUIRE(a.sign() == -1);
    REQUIRE(a.size() == 1);
    REQUIRE(a[0] == 3);
}

TEST_CASE("integer_backend - erase_first_n_words of the inline word") {
    integer_backend a(-5);
    a.erase_first_n_words(1);
    REQUIRE(a.empty());
    REQUIRE(a.sign() == 0);
    REQUIRE(!a.allocated());
}

TEST_CASE("integer_backend - erase_first_n_words past the end") {
    // erasing more words than there are used to walk _size past zero, which turned a positive value
    // into a negative one of whatever magnitude was left over
    integer_backend a {1, 2, 3};
    REQUIRE_THROWS(a.erase_first_n_words(4));

    integer_backend b {1, 2, 3};
    b.negate();
    REQUIRE_THROWS(b.erase_first_n_words(7));

    integer_backend c; // zero
    REQUIRE_THROWS(c.erase_first_n_words(1));
    c.erase_first_n_words(0); // still a no-op
    REQUIRE(c.empty());

    // erasing exactly the size is fine, and leaves zero
    integer_backend d {1, 2, 3};
    d.erase_first_n_words(3);
    REQUIRE(d.empty());
    REQUIRE(d.sign() == 0);
}

// ----- set_negative / set_zero / normalize -----

TEST_CASE("integer_backend - set_negative") {
    integer_backend a {1, 2};
    REQUIRE(a.sign() == 2);
    a.set_negative(true);
    REQUIRE(a.sign() == -2);
    REQUIRE(a.size() == 2);
    REQUIRE(a[0] == 1);
    REQUIRE(a[1] == 2);
    a.set_negative(true); // idempotent
    REQUIRE(a.sign() == -2);
    a.set_negative(false);
    REQUIRE(a.sign() == 2);
    a.set_negative(false);
    REQUIRE(a.sign() == 2);
    REQUIRE(a[1] == 2);
}

TEST_CASE("integer_backend - set_negative of zero") {
    integer_backend a;
    a.set_negative(true);
    REQUIRE(a.sign() == 0); // there is no negative zero
    REQUIRE(a.empty());
    a.set_negative(false);
    REQUIRE(a.sign() == 0);
}

TEST_CASE("integer_backend - set_negative of the inline word") {
    integer_backend a {5};
    a.set_negative(true);
    REQUIRE(a.sign() == -1);
    REQUIRE(a.size() == 1);
    REQUIRE(a[0] == 5);
    REQUIRE(!a.allocated());
}

TEST_CASE("integer_backend - set_zero of the inline word") {
    integer_backend a(-5);
    a.set_zero();
    REQUIRE(a.empty());
    REQUIRE(a.size() == 0);
    REQUIRE(a.sign() == 0);
    REQUIRE(a[0] == 0);
    REQUIRE(!a.allocated());
}

TEST_CASE("integer_backend - set_zero keeps the heap buffer") {
    integer_backend a {1, 2, 3};
    a.negate();
    const uint64_t* p = a.data();
    const int cap = a.capacity();
    a.set_zero();
    REQUIRE(a.empty());
    REQUIRE(a.sign() == 0);
    REQUIRE(a[0] == 0);
    REQUIRE(a.allocated());
    REQUIRE(a.data() == p);
    REQUIRE(a.capacity() == cap);
}

TEST_CASE("integer_backend - normalize of the inline word") {
    integer_backend a {5};
    a.normalize();
    REQUIRE(a.size() == 1);
    REQUIRE(a.sign() == 1);
    REQUIRE(a[0] == 5);

    integer_backend b {5};
    b[0] = 0;
    b.normalize();
    REQUIRE(b.empty());
    REQUIRE(b.sign() == 0);

    integer_backend c(-5);
    c[0] = 0;
    c.normalize();
    REQUIRE(c.empty());
    REQUIRE(c.sign() == 0);
}

TEST_CASE("integer_backend - normalize drops high zero words") {
    integer_backend a {7, 0, 0};
    a.normalize();
    REQUIRE(a.size() == 1);
    REQUIRE(a.sign() == 1);
    REQUIRE(a[0] == 7);
}

TEST_CASE("integer_backend - normalize of an all zero heap value") {
    integer_backend a {0, 0, 0};
    a.normalize();
    REQUIRE(a.empty());
    REQUIRE(a.size() == 0);
    REQUIRE(a.sign() == 0);
}

TEST_CASE("integer_backend - normalize keeps the sign") {
    integer_backend a {7, 0, 0};
    a.negate();
    a.normalize();
    REQUIRE(a.sign() == -1);
    REQUIRE(a.size() == 1);
    REQUIRE(a[0] == 7);

    integer_backend b {0, 0};
    b.negate();
    b.normalize();
    REQUIRE(b.sign() == 0);
    REQUIRE(b.empty());
}

TEST_CASE("integer_backend - normalize keeps interior zeros") {
    integer_backend a {0, 0, 7};
    a.normalize();
    REQUIRE(a.size() == 3);
    REQUIRE(a[0] == 0);

    integer_backend b {0, 7, 0};
    b.normalize();
    REQUIRE(b.size() == 2);
    REQUIRE(b[1] == 7);

    integer_backend c {0, 7, 0};
    c.negate();
    c.normalize();
    REQUIRE(c.sign() == -2);
}

// ----- reset_one_without_init / reset_two_without_init -----

TEST_CASE("integer_backend - reset_one_without_init") {
    integer_backend a {4, 5, 6};
    const int cap = a.capacity();
    a.reset_one_without_init();
    REQUIRE(a.size() == 1);
    REQUIRE(a.sign() == 1);
    REQUIRE(a.allocated()); // the buffer is kept
    REQUIRE(a.capacity() == cap);
    a[0] = 3;
    REQUIRE(a[0] == 3);
    REQUIRE(a.back() == 3);

    integer_backend b;
    b.reset_one_without_init();
    REQUIRE(b.size() == 1);
    REQUIRE(!b.allocated());
    b[0] = 3;
    REQUIRE(b[0] == 3);

    integer_backend c(-9);
    c.reset_one_without_init();
    REQUIRE(c.sign() == 1); // the sign is cleared
    REQUIRE(c.size() == 1);
}

TEST_CASE("integer_backend - reset_two_without_init from inline") {
    integer_backend a {5};
    a.reset_two_without_init();
    REQUIRE(a.allocated());
    REQUIRE(a.capacity() >= 2);
    REQUIRE(a.size() == 2);
    REQUIRE(a.sign() == 2);
    a[0] = 1;
    a[1] = 2;
    REQUIRE(a[0] == 1);
    REQUIRE(a[1] == 2);
    REQUIRE(a.back() == 2);
}

TEST_CASE("integer_backend - reset_two_without_init reuses the heap buffer") {
    integer_backend a {1, 2, 3, 4};
    a.negate();
    const uint64_t* p = a.data();
    const int cap = a.capacity();
    a.reset_two_without_init();
    REQUIRE(a.size() == 2);
    REQUIRE(a.sign() == 2); // the sign is cleared
    REQUIRE(a.data() == p);
    REQUIRE(a.capacity() == cap);
    a[0] = 8;
    a[1] = 9;
    REQUIRE(a[0] == 8);
    REQUIRE(a[1] == 9);
}

TEST_CASE("integer_backend - reset_two_without_init from an empty value") {
    integer_backend a;
    a.reset_two_without_init();
    REQUIRE(a.allocated());
    REQUIRE(a.size() == 2);
    a[0] = 1;
    a[1] = 0;
    a.normalize();
    REQUIRE(a.size() == 1);
}

// ----- crossing the inline/heap boundary -----

TEST_CASE("integer_backend - grow to the heap and shrink back down") {
    integer_backend a(-5);
    REQUIRE(!a.allocated());
    a.resize(4);
    REQUIRE(a.allocated());
    REQUIRE(a.sign() == -4);
    a.downsize(1);
    a.normalize();
    REQUIRE(a.sign() == -1);
    REQUIRE(a.size() == 1);
    REQUIRE(a[0] == 5);
    REQUIRE(a.allocated()); // the buffer is not released
    REQUIRE(a.capacity() >= 4);

    // and the value is still usable after the round trip
    a.insert_first_word(6);
    REQUIRE(a.sign() == -2);
    REQUIRE(a[0] == 6);
    REQUIRE(a[1] == 5);
    a.erase_first_n_words(1);
    REQUIRE(a.sign() == -1);
    REQUIRE(a[0] == 5);
}

TEST_CASE("integer_backend - insert and erase are inverses") {
    integer_backend a {1, 2, 3};
    a.negate();
    a.insert_first_n_words(2);
    a.erase_first_n_words(2);
    REQUIRE(a.sign() == -3);
    REQUIRE(a.size() == 3);
    REQUIRE(a[0] == 1);
    REQUIRE(a[1] == 2);
    REQUIRE(a[2] == 3);
}
