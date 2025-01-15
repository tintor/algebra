#include "algebra/integer_backend.h"

namespace algebra {

void integer_backend::reset(size_type size, bool initialize) {
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

void integer_backend::operator+=(word a) {
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

void integer_backend::erase_first_n_words(int n) {
    if (n > 0) {
        for (size_type i = n; i < abs(_size); i++)
            operator[](i - n) = operator[](i);
        _size += (_size >= 0) ? -n : n;
    }
}

void integer_backend::reserve(size_type capacity) {
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

void integer_backend::resize(size_type size) {
    reserve(size);
    if (_capacity)
        for (size_type i = abs(_size); i < size; i++)
            _words[i] = 0;
    _size = (_size >= 0) ? size : -size;
}

void integer_backend::insert_first_n_words(size_type n) {
    if (n == 0)
        return;

    reserve(abs(_size) + n);
    for (size_type i = abs(_size); i-- > 0; )
        operator[](i + n) = operator[](i);
    for (size_type i = 0; i < n; i++)
        operator[](i) = 0;
    _size += (_size >= 0) ? n : -n;
}

void integer_backend::insert_first_word(word a) {
    reserve(abs(_size) + 1);
    for (size_type i = abs(_size); i-- > 0; )
        operator[](i + 1) = operator[](i);
    operator[](0) = a;
    _size += (_size >= 0) ? 1 : -1;
}

}
