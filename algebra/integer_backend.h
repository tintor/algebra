#pragma once
#include <concepts>
#include <limits>
#include <cmath>

namespace algebra {

template<typename T>
auto signum(T a) { return (a > 0) - (a < 0); }

class integer_backend {
public:
    using size_type = int;
    using word = unsigned long;
    using dword = unsigned __int128;

private:
    union {
        word* _words; // least significant first
        word _single_word;
    };
    size_type _size; // negative value for negative number
    size_type _capacity;

public:
    constexpr integer_backend() : _single_word(0), _size(0), _capacity(0) { }
#if 0
    constexpr integer_backend(std::signed_integral auto a) : _single_word(static_cast<ulong>(a)), _size(signum(a)), _capacity(0) {
        // because abs(a) doesn't work for INT32_MIN  and INT64_MIN
        if (_single_word & (1lu << 63))
            _single_word = ~(_single_word - 1);
    }
#endif
    constexpr integer_backend(std::unsigned_integral auto a) : _single_word(a), _size(a > 0), _capacity(0) { }

#if 0
    constexpr integer_backend(cent a) {
        if (a > 0) {
            if (a <= std::numeric_limits<ulong>::max()) {
                _single_word = a;
                _size = 1;
                _capacity = 0;
                return;
            }
            _words = new word[2];
            _words[0] = a;
            _words[1] = a >> 64;
            _size = 2;
            _capacity = 2;
            return;
        }
        if (a < 0) {
            a = -a;
            if (a <= std::numeric_limits<ulong>::max()) {
                _single_word = a;
                _size = -1;
                _capacity = 0;
                return;
            }
            _words = new word[2];
            _words[0] = a;
            _words[1] = a >> 64;
            _size = 2;
            _capacity = 2;
            return;
        }
        _single_word = 0;
        _size = 0;
        _capacity = 0;
    }
#endif

    constexpr integer_backend(unsigned __int128 a) {
        if (a > 0) {
            if (a <= std::numeric_limits<word>::max()) {
                _single_word = a;
                _size = 1;
                _capacity = 0;
                return;
            }
            _words = new word[2];
            _words[0] = a;
            _words[1] = a >> 64;
            _size = 2;
            _capacity = 2;
            return;
        }
        _single_word = 0;
        _size = 0;
        _capacity = 0;
    }

    constexpr integer_backend(integer_backend&& o) : _words(o._words), _size(o._size), _capacity(o._capacity) {
        o._single_word = 0;
        o._size = 0;
        o._capacity = 0;
    }

    constexpr integer_backend(const integer_backend& o) : _words(o._words), _size(o._size), _capacity(o._capacity) {
        if (_capacity) {
            _words = new word[_capacity];
            for (size_type i = 0; i < _capacity; i++)
                _words[i] = o._words[i];
        }
    }

    constexpr ~integer_backend() {
        if (_capacity)
            delete[] _words;
    }

    // TODO this might not work for cent / ucent
    constexpr void operator=(std::unsigned_integral auto a) {
        if (a > 0) {
            operator[](0) = a;
            _size = 1;
            return;
        }
#if 0
        if (a < 0) {
            operator[](0) = -a;
            _size = -1;
            return;
        }
#endif
        operator[](0) = a;
        _size = 0;
    }

    constexpr void operator=(integer_backend&& o) {
        _words = o._words;
        _size = o._size;
        _capacity = o._capacity;
        o._words = 0;
        o._size = 0;
        o._capacity = 0;
    }

    void operator=(const integer_backend& o) {
        reset(o._size); // TODO unnecessary clearing inside reset()
        for (size_type i = 0; i < abs(o._size); i++)
            operator[](i) = o.operator[](i);
    }

    void reset(size_type size, bool initialize = true);
    void operator+=(word a);

    constexpr void pop_back() {
        if (_size > 0) _size--; else if (_size < 0) _size++;
    }

    constexpr void set_zero() {
        _size = 0;
        operator[](0) = 0;
    }

    size_type size() const { return abs(_size); }
    bool allocated() const { return _capacity; }
    word operator[](size_type i) const { return _capacity ? _words[i] : _single_word; }
    word& operator[](size_type i) { return _capacity ? _words[i] : _single_word; }
    word back() const { return operator[](abs(_size) - 1); }
    word& back() { return operator[](abs(_size) - 1); }

    constexpr void swap(integer_backend& o) {
        std::swap(_words, o._words);
        std::swap(_size, o._size);
        std::swap(_capacity, o._capacity);
    }

    void erase_first_n_words(size_type n);
    void reserve(size_type capacity);
    void resize(size_type size);
    void insert_first_n_words(size_type n);
    void insert_first_word(word a);

    constexpr void normalize() {
        if (_capacity) {
            if (_size > 0) {
                while (_size && _words[_size - 1] == 0)
                    _size -= 1;
            }
            if (_size < 0) {
                while (_size && _words[-_size - 1] == 0)
                    _size += 1;
            }
        } else {
            if (_single_word == 0)
                _size = 0;
        }
    }

    // part of integer backend
    void negate() { _size = -_size; }
    void set_negative(bool a) { _size = a ? -abs(_size) : abs(_size); }
    size_type sign() const { return _size; }
};

}
