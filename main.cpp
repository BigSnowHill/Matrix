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

//#include "matrix.h"
//#include "parallel_matrix.h"
#include "parallel_matrix_async.h"


struct Func{
public:
//    Func(){};
    virtual void do_func(){};
};

template<class T>
void function_clock(T f) {
    auto start = std::chrono::steady_clock::now();
    auto end = std::chrono::steady_clock::now();
    auto ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);

    start = std::chrono::steady_clock::now();

    f.do_func();

    end = std::chrono::steady_clock::now();
    ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    std::cout << ns.count() << " ns\n";
}


struct FuncEMatrix : public Func {
private:
    size_t N;
public:
    explicit FuncEMatrix(size_t n) : N(n) {};
    void do_func(){
        Matrix<double> A = Matrix<double>::E_matrix(N);
    }
};

struct FuncAddition : public Func {
private:
    Matrix<double> A, B;
public:
    explicit FuncAddition(Matrix<double> A, Matrix<double> B) : A(A), B(B) {};
    void do_func(){
        A + B;
    }
};

struct FuncMultiplication : public Func {
private:
    Matrix<double> A, B;
public:
    explicit FuncMultiplication(Matrix<double> A, Matrix<double> B) : A(A), B(B) {};
    void do_func(){
        A * B;
    }
};

struct FuncInverseMatrix : public Func{
private:
    Matrix<double> A;
public:
    explicit FuncInverseMatrix(Matrix<double> A) : A(A){};
    void do_func(){
        !A;
    }
};

int main() {
    for (size_t i = 1; i <= 1024; i*=2) {
//    for (size_t i = 1; i <= 10; ++i) {
//    size_t i = 10;
        std::cout << "Time for matrix with demension: " << i << "x" << i << std::endl;
        Matrix<double> A = Matrix<double>::E_matrix(i), B = Matrix<double>::E_matrix(i);
        FuncEMatrix f1(i); FuncAddition f2(A, B); FuncMultiplication f3(A, B); FuncInverseMatrix f4(A);
//        std::cout << "Identity matrix:"; function_clock(f1);
//        std::cout << "Matrix addition:"; function_clock(f2);
//        std::cout << "Matrix multiplication:"; function_clock(f3);
        std::cout << "Inverse matrix:"; function_clock(f4);
    }
    return 0;
}
