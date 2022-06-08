
#ifndef VECTOR_FUNCS
#define VECTOR_FUNCS

#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <stdexcept>


template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v) { 
    os << "[";
    for (int i = 0; i < v.size(); ++i) { 
        os << v[i];
        if (i != v.size() - 1)
            os << ", ";
    } 
    os << "]\n";
    return os;
} 
 
template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<std::vector<T>>& m) { 
    for (const auto& row : m) { 
        os << row;
    } 
    return os;
} 

template<typename T>
void operator+=(std::vector<T>& vec, const T& val) {
    for (T& elem : vec) {
        elem += val;
    }
}

template<typename T>
void operator+=(std::vector<std::vector<T>>& mat, const T& val) {
    for (auto& row : mat) {
        row += val;
    }
}

template<typename T>
void operator-=(std::vector<T>& vec, const T& val) {
    for (T& elem : vec) {
        elem -= val;
    }
}

template<typename T>
void operator-=(std::vector<std::vector<T>>& mat, const T& val) {
    for (auto& row : mat) {
        row -= val;
    }
}

template<typename T>
void operator-=(std::vector<T>& vec1, const std::vector<T>& vec2) {
    if (vec1.size() != vec2.size()) {
        std::ostringstream ss;
        ss << "Incompatible sizes: " << vec1.size() << " and " << vec2.size();
        throw std::invalid_argument(ss.str());
    }
    for (unsigned i=0; i<vec1.size(); i++) {
        vec1[i] -= vec2[i];
    }
}

template<typename T>
void operator-=(std::vector<std::vector<T>>& mat1, const std::vector<std::vector<T>>& mat2) {
    if (mat1.size() != mat2.size()) {
        std::ostringstream ss;
        ss << "Incompatible sizes: " << mat1.size() << " and " << mat2.size();
        throw std::invalid_argument(ss.str());
    }
    for (unsigned i=0; i<mat1.size(); i++) {
        mat1[i] -= mat2[i];
    }
}

template<typename T>
std::vector<T> subset_vector(
        const std::vector<T>& vec, 
        const std::vector<int>& cond) {
    if (vec.size() != cond.size()) {
        std::ostringstream ss;
        ss << "Incompatible sizes: " << vec.size() << " and " << cond.size();
        throw std::invalid_argument(ss.str());
    }

    std::vector<T> new_vec;
    for (unsigned i=0; i<vec.size(); i++) {
        if (cond[i] > 0) {
            new_vec.push_back(vec[i]);
        }
    }
    return new_vec;
}

template<typename T>
std::vector<std::vector<T>> select_cols(
        const std::vector<std::vector<T>>& mat, 
        const std::vector<int>& cols_idxs) {

    std::vector<std::vector<T>> new_mat;
    for (unsigned i=0; i<mat.size(); i++){
        std::vector<T> temp;
        for (const auto& j : cols_idxs) {
            temp.push_back(mat[i][j]);
        }
        new_mat.push_back(temp);
    }
    return new_mat;
}

template<typename T>
std::vector<std::vector<T>> subset_cols(
        const std::vector<std::vector<T>>& mat,
        const std::vector<int>& cond) {

    if (mat[0].size() != cond.size()) {
        std::ostringstream ss;
        ss << "Incompatible sizes: " << mat[0].size() << " and " << cond.size();
        throw std::invalid_argument(ss.str());
    }

    std::vector<std::vector<T>> new_mat;
    for (const auto& row : mat) {
        std::vector<T> temp;
        for (unsigned i=0; i < row.size(); i++) {
            if (cond[i] > 0) {
                temp.push_back(row[i]);
            }
        }
        new_mat.push_back(temp);
    }
    return new_mat;
}

template<typename T>
std::vector<std::vector<T>> subset_rows(
        const std::vector<std::vector<T>>& mat,
        const std::vector<int>& cond) {

    if (mat.size() != cond.size()) {
        std::ostringstream ss;
        ss << "Incompatible sizes: " << mat.size() << " and " << cond.size();
        throw std::invalid_argument(ss.str());
    }

    std::vector<std::vector<T>> new_mat;
    for (unsigned i=0; i<mat.size(); i++) {
        if (cond[i] > 0) {
            new_mat.push_back(mat[i]);
        }
    }
    return new_mat;
}

#endif
