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
    constexpr integer_backend(std::unsigned_integral auto a) : _single_word(a), _size((a > 0) ? 1 : 0), _capacity(0) { }
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

    constexpr void operator=(std::unsigned_integral auto a) {
        if (a > 0) {
            if constexpr (sizeof(a) == 16) {
                if (a > std::numeric_limits<uint64_t>::max()) {
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
        _words = o._words;
        _size = o._size;
        _capacity = o._capacity;
        o._words = 0;
        o._size = 0;
        o._capacity = 0;
    }

    constexpr void operator=(const integer_backend& o) {
        reset(o._size); // TODO unnecessary clearing inside reset()
        for (size_type i = 0; i < abs(o._size); i++)
            operator[](i) = o.operator[](i);
    }

    constexpr void reset(size_type size, bool initialize = true);
    constexpr void operator+=(word a);

    constexpr void pop_back() {
        if (_size > 0) _size--; else if (_size < 0) _size++;
    }

    constexpr void set_zero() {
        _size = 0;
        operator[](0) = 0;
    }

    constexpr size_type size() const { return abs(_size); }
    constexpr bool allocated() const { return _capacity; }
    constexpr word operator[](size_type i) const { return _capacity ? _words[i] : _single_word; }
    constexpr word& operator[](size_type i) { return _capacity ? _words[i] : _single_word; }
    constexpr word back() const { return operator[](abs(_size) - 1); }
    constexpr word& back() { return operator[](abs(_size) - 1); }

    constexpr void swap(integer_backend& o) {
        std::swap(_words, o._words);
        std::swap(_size, o._size);
        std::swap(_capacity, o._capacity);
    }

    constexpr void erase_first_n_words(size_type n);
    constexpr void reserve(size_type capacity);
    constexpr void resize(size_type size);
    constexpr void insert_first_n_words(size_type n);
    constexpr void insert_first_word(word a);

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
    constexpr void set_negative(bool a) { _size = a ? -abs(_size) : abs(_size); }
    constexpr size_type sign() const { return _size; }
};

constexpr void integer_backend::reset(size_type size, bool initialize) {
    if (_capacity) {
        // heap
        if (abs(size) > _capacity) {
            // heap -> heap x2
            _capacity = std::max(_capacity * 2, abs(size));
            delete[] _words;
            _words = new word[_capacity];
        } else if (abs(size) <= 1) {
            // heap -> SBO
            delete[] _words;
            _single_word = 0;
            _capacity = 0;
        }
    } else {
        // SBO
        if (abs(size) <= 1) {
            _single_word = 0;
        } else {
            // SBO -> heap
            _capacity = abs(size);
            _words = new word[_capacity];
        }
    }
    _size = size;

    if (initialize && _capacity)
        for (size_type i = 0; i < abs(size); i++)
            _words[i] = 0;
}

constexpr void integer_backend::operator+=(word a) {
    if (_capacity) {
        if (abs(_size) == _capacity) {
            // heap -> heap x2
            word* w = new word[_capacity * 2];
            for (word i = 0; i < _capacity; i++)
                w[i] = _words[i];
            for (word i = 0; i < _capacity; i++)
                w[_capacity + i] = 0;
            delete[] _words;
            _words = w;
            _capacity *= 2;
        }

        if (_size < _capacity) {
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
            word* w = new word[2];
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
        for (size_type i = n; i < abs(_size); i++)
            operator[](i - n) = operator[](i);
        _size += (_size >= 0) ? -n : n;
    }
}

constexpr void integer_backend::reserve(size_type capacity) {
    if (capacity <= _capacity || capacity == 1)
        return;
    capacity = std::max(_capacity * 2, capacity);
    word* words = new word[capacity];
    for (size_type i = 0; i < abs(_size); i++)
        words[i] = operator[](i);
    if (_capacity)
        delete[] _words;
    _words = words;
    _capacity = capacity;
}

constexpr void integer_backend::resize(size_type size) {
    reserve(size);
    if (_capacity)
        for (size_type i = abs(_size); i < size; i++)
            _words[i] = 0;
    _size = (_size >= 0) ? size : -size;
}

constexpr void integer_backend::insert_first_n_words(size_type n) {
    if (n == 0)
        return;

    reserve(abs(_size) + n);
    for (size_type i = abs(_size); i-- > 0; )
        operator[](i + n) = operator[](i);
    for (size_type i = 0; i < n; i++)
        operator[](i) = 0;
    _size += (_size >= 0) ? n : -n;
}

constexpr void integer_backend::insert_first_word(word a) {
    reserve(abs(_size) + 1);
    for (size_type i = abs(_size); i-- > 0; )
        operator[](i + 1) = operator[](i);
    operator[](0) = a;
    _size += (_size >= 0) ? 1 : -1;
}

}
