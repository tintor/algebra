#pragma once
#include "algebra/types.h"

namespace algebra {

template<std::floating_point T> auto signum(T a) { return (a > 0) - (a < 0); }
template<std_int T> auto signum(T a) { return (a > 0) - (a < 0); }

class integer_backend {
private:
    union {
        uint64_t* _words; // least significant first
        uint64_t _single_word;
    };
    int _size; // negative value for negative number
    int _capacity;

public:
    constexpr integer_backend(std::initializer_list<uint64_t> a) : _words(0), _size(0), _capacity(0) {
        resize(a.size());
        for (int i = 0; i < a.size(); i++)
            operator[](i) = a.begin()[i];
    }
    constexpr integer_backend() : _single_word(0), _size(0), _capacity(0) { }
    constexpr integer_backend(std_int auto a) {
        const auto au = abs_unsigned(a);

        static_assert(sizeof(a) == 16 || sizeof(a) <= 8);
        if (a == 0) {
            _single_word = 0;
            _size = 0;
            _capacity = 0;
            return;
        }
        if (au <= UINT64_MAX) {
            _single_word = au;
            _size = (a < 0) ? -1 : 1;
            _capacity = 0;
            return;
        }
        if constexpr (sizeof(a) == 16) {
            _words = new uint64_t[2];
            _words[0] = au;
            _words[1] = au >> 64;
            _size = (a < 0) ? -2 : 2;
            _capacity = 2;
        }
    }

    constexpr integer_backend(integer_backend&& o) : _words(o._words), _size(o._size), _capacity(o._capacity) {
        o._single_word = 0;
        o._size = 0;
        o._capacity = 0;
    }

    constexpr integer_backend(const integer_backend& o) {
        if (o.size() >= 2) {
            _capacity = o.size();
            _words = new uint64_t[_capacity];
            for (int i = 0; i < _capacity; i++)
                _words[i] = o._words[i];
        } else {
            _capacity = 0;
            _single_word = o[0];
        }
        _size = o._size;
    }

    constexpr ~integer_backend() {
        if (_capacity)
            delete[] _words;
    }

    constexpr void operator=(std_int auto a) {
        const auto au = abs_unsigned(a);
        static_assert(sizeof(a) == 16 || sizeof(a) <= 8);

        if (a == 0) {
            set_zero();
            return;
        }
        if (au <= UINT64_MAX) {
            _size = (a < 0) ? -1 : 1;
            operator[](0) = au;
            return;
        }
        if constexpr (sizeof(au) == 16) {
            reset_two_without_init();
            _size = (a < 0) ? -2 : 2;
            operator[](0) = au;
            operator[](1) = au >> 64;
        }
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
        for (int i = 0; i < size(); i++)
            operator[](i) = o.operator[](i);
    }

    constexpr void reset_one_without_init() { _size = 1; }
    constexpr void reset_two_without_init() {
        if (!_capacity) {
            _capacity = 2;
            _words = new uint64_t[2];
        }
        _size = 2;
    }
    constexpr void reset(int size, bool initialize = true);
    constexpr void push_back(uint64_t a);
    constexpr void push_back(uint64_t a, uint64_t b);

    constexpr void pop_back() {
        if (_size > 0) _size--; else if (_size < 0) _size++;
    }

    constexpr void set_zero() {
        _size = 0;
        operator[](0) = 0;
    }

    constexpr int capacity() const { return std::max(1, _capacity); }
    constexpr int size() const { return std::abs(_size); }
    constexpr bool empty() const { return _size == 0; }
    constexpr bool allocated() const { return _capacity; }
    constexpr const uint64_t* data() const { return _capacity ? _words : &_single_word; }
    constexpr uint64_t* data() { return _capacity ? _words : &_single_word; }
    constexpr uint64_t operator[](int i) const { return _capacity ? _words[i] : _single_word; }
    constexpr uint64_t& operator[](int i) { return _capacity ? _words[i] : _single_word; }
    constexpr uint64_t back() const { return operator[](size() - 1); }
    constexpr uint64_t& back() { return operator[](size() - 1); }

    constexpr void swap(integer_backend& o) {
        std::swap(_words, o._words);
        std::swap(_size, o._size);
        std::swap(_capacity, o._capacity);
    }

    constexpr void erase_first_n_words(int n);
    constexpr void reserve(int capacity);
    constexpr void reserve_and_set_zero(int capacity);
    constexpr void reserve_bits(size_t bits) { return reserve((bits + 63) / 64); }
    constexpr void resize(int size);
    constexpr void downsize(int size);

    constexpr void insert_first_n_words(int n);
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
    constexpr int sign() const { return _size; }

    static constexpr uint64_t hash_fn_64bit(uint64_t k) {
        k ^= k >> 33;
        k *= 0xff51afd7ed558ccdllu;
        return k;
    }
};

constexpr void integer_backend::reset(int size, bool initialize) {
    const int abs_size = std::abs(size);
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
        for (int i = 0; i < abs_size; i++)
            _words[i] = 0;
}

constexpr void integer_backend::push_back(uint64_t a) {
    if (_capacity) {
        if (size() == _capacity) {
            // heap -> heap x2
            auto w = new uint64_t[_capacity * 2];
            for (int i = 0; i < _capacity; i++)
                w[i] = _words[i];
            for (int i = 0; i < _capacity; i++)
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

constexpr void integer_backend::push_back(uint64_t a, uint64_t b) {
    // TODO optimize
    push_back(a);
    push_back(b);
}

constexpr void integer_backend::erase_first_n_words(int n) {
    if (n > 0) {
        for (int i = n; i < size(); i++)
            operator[](i - n) = operator[](i);
        _size += (_size >= 0) ? -n : n;
    }
}

constexpr void integer_backend::reserve_and_set_zero(int capacity) {
    _size = 0;
    if (capacity <= _capacity || capacity == 1) {
        operator[](0) = 0;
        return;
    }
    if (_capacity)
        delete[] _words;
    _capacity = std::max(_capacity * 2, capacity);
    _words = new uint64_t[_capacity];
}

constexpr void integer_backend::reserve(int capacity) {
    if (capacity <= _capacity || capacity == 1)
        return;
    capacity = std::max(_capacity * 2, capacity);
    auto words = new uint64_t[capacity];
    for (int i = 0; i < size(); i++)
        words[i] = operator[](i);
    if (_capacity)
        delete[] _words;
    _words = words;
    _capacity = capacity;
}

constexpr void integer_backend::resize(int size) {
    reserve(size);
    if (_capacity)
        for (int i = this->size(); i < size; i++)
            _words[i] = 0;
    _size = (_size >= 0) ? size : -size;
}

constexpr void integer_backend::downsize(int size) {
    _size = (_size >= 0) ? size : -size;
}

constexpr void integer_backend::insert_first_n_words(int n) {
    if (n == 0)
        return;

    reserve(size() + n);
    for (int i = size(); i-- > 0; )
        operator[](i + n) = operator[](i);
    for (int i = 0; i < n; i++)
        operator[](i) = 0;
    _size += (_size >= 0) ? n : -n;
}

constexpr void integer_backend::insert_first_word(uint64_t a) {
    reserve(size() + 1);
    for (int i = size(); i-- > 0; )
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
