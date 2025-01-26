#pragma once
#include <concepts>
#include <limits>
#include <cmath>
#include <type_traits>

namespace algebra {

template<std::integral T> constexpr std::make_unsigned_t<T> make_unsigned(T a) { return a; }

using int128_t = __int128;
using uint128_t = unsigned __int128;

template<typename T>
auto signum(T a) { return (a > 0) - (a < 0); }

class integer_backend {
public:
    using size_type = int;

private:
    union {
        uint64_t* _words; // least significant first
        uint64_t _single_word;
    };
    size_type _size; // negative value for negative number
    size_type _capacity;

public:
    constexpr integer_backend() : _single_word(0), _size(0), _capacity(0) { }
    constexpr integer_backend(std::unsigned_integral auto a) : _single_word(a), _size((a > 0) ? 1 : 0), _capacity(0) { }
    constexpr integer_backend(uint128_t a) {
        if (a > 0) {
            if (a <= UINT64_MAX) {
                _single_word = a;
                _size = 1;
                _capacity = 0;
                return;
            }
            _words = new uint64_t[2];
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
            _words = new uint64_t[_capacity];
            for (size_type i = 0; i < _capacity; i++)
                _words[i] = o._words[i];
        }
    }

    constexpr ~integer_backend() {
        if (_capacity)
            delete[] _words;
    }

    constexpr void operator=(std::unsigned_integral auto a) {
        if (a > 0) {
            if constexpr (sizeof(a) == 16) {
                if (a > UINT64_MAX) {
                    operator[](0) = a;
                    operator[](1) = a >> 64;
                    _size = 2;
                } else {
                    operator[](0) = a;
                    _size = 1;
                }
                return;
            }
            operator[](0) = a;
            _size = 1;
            return;
        }
        operator[](0) = 0;
        _size = 0;
    }

    constexpr void operator=(integer_backend&& o) {
        if (this == &o)
            return;
        _words = o._words;
        _size = o._size;
        _capacity = o._capacity;
        o._words = 0;
        o._size = 0;
        o._capacity = 0;
    }

    constexpr void operator=(const integer_backend& o) {
        if (this == &o)
            return;
        reset(o._size, /*initialize*/false);
        for (size_type i = 0; i < size(); i++)
            operator[](i) = o.operator[](i);
    }

    constexpr void reset_one_without_init() { _size = 1; }
    constexpr void reset(size_type size, bool initialize = true);
    constexpr void operator+=(uint64_t a);

    constexpr void pop_back() {
        if (_size > 0) _size--; else if (_size < 0) _size++;
    }

    constexpr void set_zero() {
        _size = 0;
        operator[](0) = 0;
    }

    constexpr size_type size() const { return std::abs(_size); }
    constexpr bool allocated() const { return _capacity; }
    constexpr uint64_t operator[](size_type i) const { return _capacity ? _words[i] : _single_word; }
    constexpr uint64_t& operator[](size_type i) { return _capacity ? _words[i] : _single_word; }
    constexpr uint64_t back() const { return operator[](size() - 1); }
    constexpr uint64_t& back() { return operator[](size() - 1); }

    constexpr void swap(integer_backend& o) {
        std::swap(_words, o._words);
        std::swap(_size, o._size);
        std::swap(_capacity, o._capacity);
    }

    constexpr void erase_first_n_words(size_type n);
    constexpr void reserve(size_type capacity);
    constexpr void resize(size_type size);

    constexpr void insert_first_n_words(size_type n);
    constexpr void insert_first_word(uint64_t a);

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

    // not accessible to class natural
    constexpr void negate() { _size = -_size; }
    constexpr void set_negative(bool a) { _size = a ? -size() : size(); }
    constexpr size_type sign() const { return _size; }

    static constexpr uint64_t hash_fn_64bit(uint64_t k) {
        k ^= k >> 33;
        k *= 0xff51afd7ed558ccdllu;
        return k;
    }
};

constexpr void integer_backend::reset(size_type size, bool initialize) {
    const size_type abs_size = std::abs(size);
    if (_capacity) {
        // heap
        if (abs_size > _capacity) {
            // heap -> heap x2
            _capacity = std::max(_capacity * 2, abs_size);
            delete[] _words;
            _words = new uint64_t[_capacity];
        } else if (abs_size <= 1) {
            // heap -> SBO
            delete[] _words;
            _single_word = 0;
            _capacity = 0;
        }
    } else {
        // SBO
        if (abs_size <= 1) {
            _single_word = 0;
        } else {
            // SBO -> heap
            _capacity = abs_size;
            _words = new uint64_t[_capacity];
        }
    }
    _size = size;

    if (initialize && _capacity)
        for (size_type i = 0; i < abs_size; i++)
            _words[i] = 0;
}

constexpr void integer_backend::operator+=(uint64_t a) {
    if (_capacity) {
        if (size() == _capacity) {
            // heap -> heap x2
            auto w = new uint64_t[_capacity * 2];
            for (size_type i = 0; i < _capacity; i++)
                w[i] = _words[i];
            for (size_type i = 0; i < _capacity; i++)
                w[_capacity + i] = 0;
            delete[] _words;
            _words = w;
            _capacity *= 2;
        }

        if (_size >= 0) {
            _words[_size] = a;
            _size++;
        } else {
            _words[-_size] = a;
            _size--;
        }
    } else {
        // SBO
        if (_size == 0) {
            _size = 1;
            _single_word = a;
        } else {
            // SBO -> heap
            auto w = new uint64_t[2];
            _capacity = 2;
            w[0] = _single_word;
            _words = w;
            w[1] = a;
            _size *= 2;
        }
    }
}

constexpr void integer_backend::erase_first_n_words(int n) {
    if (n > 0) {
        for (size_type i = n; i < size(); i++)
            operator[](i - n) = operator[](i);
        _size += (_size >= 0) ? -n : n;
    }
}

constexpr void integer_backend::reserve(size_type capacity) {
    if (capacity <= _capacity || capacity == 1)
        return;
    capacity = std::max(_capacity * 2, capacity);
    auto words = new uint64_t[capacity];
    for (size_type i = 0; i < size(); i++)
        words[i] = operator[](i);
    if (_capacity)
        delete[] _words;
    _words = words;
    _capacity = capacity;
}

constexpr void integer_backend::resize(size_type size) {
    reserve(size);
    if (_capacity)
        for (size_type i = this->size(); i < size; i++)
            _words[i] = 0;
    _size = (_size >= 0) ? size : -size;
}

constexpr void integer_backend::insert_first_n_words(size_type n) {
    if (n == 0)
        return;

    reserve(size() + n);
    for (size_type i = size(); i-- > 0; )
        operator[](i + n) = operator[](i);
    for (size_type i = 0; i < n; i++)
        operator[](i) = 0;
    _size += (_size >= 0) ? n : -n;
}

constexpr void integer_backend::insert_first_word(uint64_t a) {
    reserve(size() + 1);
    for (size_type i = size(); i-- > 0; )
        operator[](i + 1) = operator[](i);
    operator[](0) = a;
    _size += (_size >= 0) ? 1 : -1;
}

}

template<>
struct std::hash<algebra::integer_backend> {
    constexpr size_t operator()(const algebra::integer_backend& a) const {
        uint64_t seed = algebra::integer_backend::hash_fn_64bit(a.sign());
        for (int i = 0; i < a.size(); i++)
            seed = algebra::integer_backend::hash_fn_64bit(seed ^ a[i]);
        return seed;
    }
};
