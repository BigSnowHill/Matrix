#ifndef LAB_13_PARALLEL_MATRIX_H
#define LAB_13_PARALLEL_MATRIX_H

#include <iostream>
#include <vector>
#include <queue>
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <thread>
#include <chrono>
#include <functional>
#include <future>
//#include <pthread>

const size_t THREADS_COUNT = 2;


template<typename T>
class Matrix {
public:
private:
    size_t N, M;
    std::vector <std::vector<T>> matrix;
public:
    // Cunstractor if all elemenets are known
    explicit Matrix(size_t n, size_t m, std::vector <std::vector<T>> new_matrix) : N(n), M(m), matrix(new_matrix) {};

    // Cunstractor if only size is known, take elemenets from console
    Matrix(size_t n, size_t m) : N(n), M(m) {
        for (size_t i = 0; i < N; ++i) {
            std::vector <T> elem_vector;
            for (size_t j = 0; j < M; ++j) {
                T element;
                std::cin >> element;
                elem_vector.push_back(element);
            }
            matrix.push_back(elem_vector);
        }
    }

    // Cunstractor take size and elements from console
    Matrix() {
        std::cin >> N >> M;
        for (size_t i = 0; i < N; ++i) {
            std::vector <T> elem_vector;
            for (size_t j = 0; j < M; ++j) {
                T element;
                std::cin >> element;
                elem_vector.push_back(element);
            }
            matrix.push_back(elem_vector);
        }
    }

    // Copy Constractor
//    Matrix(Matrix<T> &A) : N(A.get_lines_number()), M(A.get_colomns_number()) {
//        std::vector <T> line(N, 0);
//        std::vector <std::vector<T>> lines(M, line);
//        for (int i = 0; i < N; ++i) {
//            for (int j = 0; j < M; ++j) {
//                lines[i][j] = A[i][j];
//            }
//        }
//        matrix = lines;
//    }

    // File Constractor
    Matrix(const std::string &filename) {
        std::ifstream file_in(filename);
        file_in >> N >> M;
        for (size_t i = 0; i < N; ++i) {
            std::vector <T> elem_vector;
            for (size_t j = 0; j < M; ++j) {
                T element;
                file_in >> element;
                elem_vector.push_back(element);
            }
            matrix.push_back(elem_vector);
        }
    }

    size_t get_lines_number() {
        return N;
    }

    size_t get_colomns_number() {
        return M;
    }

    std::vector <T> &operator[](size_t i) {
        return matrix[i];
    }

    void operator*=(T number) {
        size_t th_n = THREADS_COUNT;
        size_t el_n = (N * M) / th_n;
        Matrix<T> Res = *this;
        std::thread th[th_n];
        std::function<void(size_t)> func;
        func = [&](size_t num) {
            for (size_t j = (el_n * num); j < el_n * (num + 1); ++j) {
                Res[j / M][j % M] *= number;
            }
            if (num + 1 == th_n && (N * M) % th_n != 0) {
                for (size_t i = el_n * (num + 1); i < el_n * (num + 1) + ((N * M) % th_n); ++i) {
                    Res[i / M][i % M] *= number;
                }
            }
        };

        for (size_t i = 0; i < th_n; ++i) {
            th[i] = std::thread(func, i);
        }

        for (size_t i = 0; i < th_n; ++i) {
            th[i].join();
        }
    }

    static Matrix<T> zero_matrix(size_t n, size_t m) {
        std::vector <T> B(m, 0);
        std::vector <std::vector<T>> A(n, B);
        Matrix<T>* Zereos = new Matrix(n, m, A);
        return *Zereos;
    }


    static Matrix<T> E_matrix(size_t n) {
        size_t th_n = THREADS_COUNT;
        size_t el_n = n / th_n;
        Matrix<T> M = Matrix<T>::zero_matrix(n, n);
        std::thread th[th_n];
        std::function<void(size_t)> func;
        func = [&](size_t num) {
            for (size_t i = num * el_n; i < (num + 1) * el_n; ++i) {
                M[i][i] = 1;
            }
            if (num + 1 == th_n && n % th_n != 0) {
                for (size_t j = (num + 1) * el_n; j < n; ++j) {
                    M[j][j] = 1;
                }
            }
        };

        for (size_t i = 0; i < th_n; ++i) {
            th[i] = std::thread(func, i);
        }

        for (size_t i = 0; i < th_n; ++i) {
            th[i].join();
        }
        return M;
    }
};

template<class T>
T determinant(Matrix<T>);

template<class T>
Matrix<T> Transporation(Matrix<T> A) {
    if (A.get_lines_number() != A.get_colomns_number()) {
        std::cerr << "Impossible to make operation for non square matrix.\n";
        exit(EXIT_FAILURE);
    }
    std::vector <T> line(A.get_colomns_number(), 0);
    std::vector <std::vector<T>> lines(A.get_lines_number(), line);
    for (int i = 0; i < A.get_lines_number(); ++i) {
        for (int j = 0; j < A.get_colomns_number(); ++j) {
            lines[j][i] = A[i][j];
        }
    }
    Matrix<T> *B = new Matrix<T>(A.get_lines_number(), A.get_colomns_number(), lines);
    return *B;
}

template<class T>
T Alg_Ad(Matrix<T> A, size_t x, size_t y) {
    if (A.get_lines_number() != A.get_colomns_number()) {
        std::cerr << "Impossible to find determinante for not square matrix.\n";
        exit(EXIT_FAILURE);
    }
    if (A.get_lines_number() == 1) {
        return A[0][0];
    }
    std::vector <std::vector<T>> B_v;
    for (int i = 0; i < A.get_lines_number(); ++i) {
        if (i == x) continue;
        std::vector <T> line;
        for (int j = 0; j < A.get_colomns_number(); ++j) {
            if (j == y) continue;
            line.push_back(A[i][j]);
        }
        B_v.push_back(line);
    }
    Matrix<T> B(A.get_lines_number() - 1, A.get_colomns_number() - 1, B_v);
    return pow(-1, x + y) * determinant<T>(B);
}

template<class T>
T determinant(Matrix<T> A) {
    if (A.get_colomns_number() != A.get_lines_number()) {
        std::cerr << "Impossible to find determinante for not square matrix.\n";
        exit(EXIT_FAILURE);
    }
    T det = 0;
    size_t n = A.get_lines_number();
    size_t th_n = THREADS_COUNT;
    size_t el_n = (1 * n) / th_n;
    std::thread th[th_n];
    std::function<void(size_t)> func;
    func = [&](size_t f_th_n) {
        for (size_t i = (el_n * f_th_n); i < el_n * (f_th_n + 1); ++i) {
            det += A[0][i] * Alg_Ad<T>(A, 0, i);
        }
        if (f_th_n + 1 == th_n && (n * n) % th_n != 0) {
            for (size_t j = el_n * (f_th_n + 1); j < el_n * (f_th_n + 1) + ((n * n) % th_n); ++j) {
                det += A[0][j] * Alg_Ad(A, 0, j);
            }
        }
    };
    for (size_t i = 0; i < th_n; ++i) {
        th[i] = std::thread(func, i);
    }
    for (size_t i = 0; i < th_n; ++i) {
        th[i].join();
    }
    return det;
}

template<class T>
Matrix<double> &operator!(Matrix<T> A) {
    if (A.get_colomns_number() != A.get_lines_number()) {
        std::cerr << "Impossible to find determinante for not square matrix.\n";
        exit(EXIT_FAILURE);
    }
    size_t n = A.get_lines_number();
    size_t th_n = THREADS_COUNT;
    size_t el_n = (n * n) / th_n;
    std::thread th[th_n];
    Matrix<double> B = Matrix<T>::zero_matrix(n, n);
    std::function<void(size_t)> func;
    func = [&](size_t f_th_n) {
        for (size_t i = (el_n * f_th_n); i < el_n * (f_th_n + 1); ++i) {
            B[i / n][i % n] = Alg_Ad(A, i / n, i % n);
        }
        if (f_th_n + 1 == th_n && (n * n) % th_n != 0) {
            for (size_t j = el_n * (f_th_n + 1); j < f_th_n * (f_th_n + 1) + ((n * n) % f_th_n); ++j) {
                B[j / n][j % n] = Alg_Ad(A, j / n, j % n);
            }
        }
    };

    for (size_t i = 0; i < th_n; ++i) {
        th[i] = std::thread(func, i);
    }

    for (size_t i = 0; i < th_n; ++i) {
        th[i].join();
    }
    B *= (1 / determinant(A));
    Matrix<T> *C = new Matrix<T>(Transporation(B));
    return *(C);
}


template<class T>
std::ostream &operator<<(std::ostream &out, Matrix<T> A) {
    for (size_t i = 0; i < A.get_lines_number(); ++i) {
        for (size_t j = 0; j < A.get_colomns_number(); ++j) {
            out << A[i][j] << " ";
        }
        out << "\n";
    }
    return out;
}


template<class T>
Matrix<T> &operator*(Matrix<T> A, Matrix<T> B) {
    if (A.get_colomns_number() != B.get_lines_number()) {
        std::cerr << "Error: impossible to multiply matrices.\n";
        exit(EXIT_FAILURE);
    }
    size_t n = A.get_lines_number(), m = A.get_colomns_number();
    size_t th_n = THREADS_COUNT;
    size_t el_n = (n * m) / th_n;
    Matrix<T> Res = Matrix<T>::zero_matrix(m, n);
    std::thread th[th_n];
    std::function < T(Matrix<T> &, Matrix<T>&, size_t, size_t)> el_func;
    std::function<void(Matrix<T> &, Matrix<T> &, Matrix<T>&, size_t)> func;

    el_func = [](Matrix<T> &N, Matrix<T> &M, size_t line, size_t column) {
        T sum = 0;
        for (size_t i = 0; i < column; ++i) {
            sum += N[line][i] * M[i][column];
        }
        return sum;
    };

    func = [&](Matrix<T> &Res, Matrix<T> &C, Matrix<T> &F, size_t f_th_n) {
        for (size_t j = (el_n * th_n); j < el_n * (th_n + 1); ++j) {
            Res[j / m][j % m] = el_func(C, F, j / m, j % m);
        }
        if (f_th_n + 1 == th_n && (n * m) % th_n != 0) {
            for (size_t i = el_n * (f_th_n + 1); i < el_n * (f_th_n + 1) + ((n * m) % th_n); ++i) {
                Res[i / m][i % m] = el_func(C, F, i / m, i % m);
            }
        }
    };
    for (size_t i = 0; i < th_n; ++i) {
        th[i] = std::thread(func, std::ref(Res), std::ref(A), std::ref(B), i);
    }
    for (size_t i = 0; i < th_n; ++i) {
        th[i].join();
    }
    Matrix<T>* Res1 = new Matrix<T>(Res);
    return *Res1;
}


template<class T>
Matrix<T> &operator*(Matrix<T> A, T c) {
    size_t n = A.get_lines_number(), m = A.get_colomns_number();
    size_t th_n = THREADS_COUNT;
    size_t el_n = (n * m) / th_n;
    Matrix<T> Res = A;
    std::thread th[th_n];
    std::function<void(size_t)> func;
    func = [&](size_t num) {
        for (size_t j = (el_n * num); j < el_n * (num + 1); ++j) {
            Res[j / m][j % m] *= c;
        }
        if (num + 1 == th_n && (n * m) % th_n != 0) {
            for (size_t i = el_n * (num + 1); i < el_n * (num + 1) + ((n * m) % th_n); ++i) {
                Res[i / m][i % m] *= c;
            }
        }
    };
    for (size_t i = 0; i < th_n; ++i) {
        th[i] = std::thread(func, i);
    }
    for (size_t i = 0; i < th_n; ++i) {
        th[i].join();
    }
    return Res;
}

template<class T>
Matrix<T> &operator*(T c, Matrix<T> A){
    return A * c;
}


template<class T>
Matrix<T> operator+(Matrix<T> &M, Matrix<T> &N) {
    if (M.get_lines_number() != N.get_lines_number() || M.get_colomns_number() != M.get_lines_number()) {
        std::cerr << "Error: impossible to add matrices.";
        exit(EXIT_FAILURE);
    }
    size_t n = N.get_lines_number(), m = N.get_colomns_number();
    size_t th_n = THREADS_COUNT;
    size_t el_n = (n * m) / th_n;
    Matrix<T> Res = N;
    std::thread th[th_n];
    std::function<void(size_t)> func;
    func = [&](size_t f_th_n) {
        for (size_t i = (el_n * f_th_n); i < el_n * (f_th_n + 1); ++i) {
            Res[i / m][i % m] += M[i / m][i % m];
        }
        if (f_th_n + 1 == th_n && (n * m) % th_n != 0) {
            for (size_t j = el_n * (f_th_n + 1); j < el_n * (f_th_n + 1) + ((n * m) % th_n); ++j) {
                Res[j / m][j % m] += M[j / m][j % m];
            }
        }
    };

    for (size_t i = 0; i < th_n; ++i) {
        th[i] = std::thread(func, i);
    }

    for (size_t i = 0; i < th_n; ++i) {
        th[i].join();
    }
    return Res;
}


template<class T>
Matrix<T> operator-(Matrix<T> &N, Matrix<T> &M) {
    if (N.get_lines_number() != M.get_lines_number || N.get_colomns_number() != M.get_colomns_number()) {
        std::cerr << "Error: impossible to subtract matrices.";
        exit(EXIT_FAILURE);
    }

    size_t n = N.get_lines_number(), m = N.get_colomns_number();
    size_t th_n = THREADS_COUNT;
    size_t el_n = (n * m) / th_n;
    Matrix<T> Res = N;
    std::thread th[th_n];
    std::function<void(size_t)> func;

    func = [&](size_t f_th_n) {
        for (size_t i = (el_n * f_th_n); i < el_n * (f_th_n + 1); ++i) {
            Res[i / m][i % m] -= M[i / m][i % m];
        }
        if (f_th_n + 1 == th_n && (n * m) % th_n != 0) {
            for (size_t j = el_n * (f_th_n + 1); j < el_n * (f_th_n + 1) + ((n * m) % th_n); ++j) {
                Res[j / m][j % m] -= M[j / m][j % m];
            }
        }
    };

    for (size_t i = 0; i < th_n; ++i) {
        th[i] = std::thread(func, i);
    }
    for (size_t i = 0; i < th_n; ++i) {
        th[i].join();
    }
    return Res;
}


template<class T>
bool operator==(Matrix<T> A, Matrix<T> B) {
    if (A.get_colomns_number() != B.get_colomns_number() || A.get_lines_number() != B.get_lines_number()) {
        return false;
    }
    for (int i = 0; i < A.get_lines_number(); ++i) {
        for (int j = 0; j < B.get_colomns_number(); ++j) {
            if (A[i][j] != B[i][j]) {
                return false;
            }
        }
    }
    return true;
}


template<class T>
bool operator!=(Matrix<T> A, Matrix<T> B) {
    if (A.get_colomns_number() != B.get_colomns_number() || A.get_lines_number() != B.get_lines_number()) {
        return true;
    }
    for (int i = 0; i < A.get_lines_number(); ++i) {
        for (int j = 0; j < B.get_colomns_number(); ++j) {
            if (A[i][j] == B[i][j]) {
                return false;
            }
        }
    }
    return true;
}


template<class T>
bool operator==(Matrix<T> &A, T c) {
    size_t n = A.get_lines_number(), m = A.get_colomns_number();
    if (c == 0) {
        Matrix<T> B = Matrix<T>::zero_matrix(n, m);
        return A == B;
    }
    if (n != m) {
        std::cerr << "Impossible to compare not square matrix and scalar.";
        exit(EXIT_FAILURE);
    }
    Matrix<T> B = Matrix<T>::E_matrix(n);
    B *= c;
    return A == B;
}


template<class T>
bool operator==(T c, Matrix<T> &A) {
    size_t n = A.get_lines_number(), m = A.get_colomns_number();
    if (c == 0) {
        Matrix<T> B = Matrix<T>::zero_matrix(n, m);
        return A == B;
    }
    if (n != m) {
        std::cerr << "Impossible to compare not square matrix and scalar.";
        exit(EXIT_FAILURE);
    }
    Matrix<T> B = Matrix<T>::E_matrix(n);
    B *= c;
    return A == B;
}


#endif //LAB_13_PARALLEL_MATRIX_H
