#ifndef MYMATRIX_H
#define MYMATRIX_H

#include "myvector.hpp"
#include "util.hpp"

template <typename T>
class mymatrix;

template <typename T>
std::ostream &operator<<(std::ostream &os, const mymatrix<T> &vec);

template <typename T>
mymatrix<T> mydot(const mymatrix<T> &mat1, const mymatrix<T> &mat2);

template <typename T>
mymatrix<T> myouter(const mymatrix<T> &mat1, const mymatrix<T> &mat2);

template <typename T>
double mynorm(const mymatrix<T> &mat, const char *ord = "L2");

template <class T>
class mymatrix
{
private:
    int n_, m_;
    T **v;

public:
    mymatrix();
    mymatrix(int n, int m);
    mymatrix(int n, int m, const T &a);
    mymatrix(int n, int m, const std::vector<std::vector<T>> a);
    mymatrix(const mymatrix &mat);

    inline mymatrix &operator=(const mymatrix &mat);
    inline T *operator[](const int i);
    inline const T *operator[](const int i) const;
    inline mymatrix operator+(const mymatrix &mat) const;
    inline mymatrix operator-(const mymatrix &mat) const;
    inline mymatrix &operator+=(const mymatrix &mat);
    inline mymatrix &operator-=(const mymatrix &mat);
    friend std::ostream &operator<<<T>(std::ostream &os, const mymatrix<T> &vec);

    inline int nrows() const;
    inline int ncols() const;
    void resize(int newn, int newm);
    void assign(int newn, int newm, const T &a);
    void identity();
    friend mymatrix<T> mydot<T>(const mymatrix<T> &mat1, const mymatrix<T> &mat2);
    friend mymatrix<T> myouter<T>(const mymatrix<T> &mat1, const mymatrix<T> &mat2);
    friend double mynorm<T>(const mymatrix<T> &mat, const char *ord);

    ~mymatrix();
};

template <class T>
mymatrix<T>::mymatrix() : n_(0), m_(0), v(nullptr) {}

template <class T>
mymatrix<T>::mymatrix(int n, int m) : n_(n), m_(m), v(n > 0 ? new T *[n] : nullptr)
{
    int size = n * m;
    if (v)
        v[0] = size > 0 ? new T[size] : nullptr;
    for (int i = 1; i < n; ++i)
        v[i] = v[i - 1] + m;
}

template <class T>
mymatrix<T>::mymatrix(int n, int m, const T &a) : n_(n), m_(m), v(n > 0 ? new T *[n] : nullptr)
{
    int size = n * m;
    if (v)
        v[0] = size > 0 ? new T[size] : nullptr;
    for (int i = 1; i < n; ++i)
        v[i] = v[i - 1] + m;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < m; ++j)
            v[i][j] = a;
}

template <class T>
mymatrix<T>::mymatrix(int n, int m, const std::vector<std::vector<T>> a) : n_(n), m_(m), v(n > 0 ? new T *[n] : nullptr)
{
    int size = n * m;
    if (v)
        v[0] = size > 0 ? new T[size] : nullptr;
    for (int i = 1; i < n; ++i)
        v[i] = v[i - 1] + m;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < m; ++j)
            v[i][j] = a[i][j];
}

template <class T>
mymatrix<T>::mymatrix(const mymatrix &mat) : n_(mat.n_), m_(mat.m_), v(n_ > 0 ? new T *[n_] : nullptr)
{
    int size = n_ * m_;
    if (v)
        v[0] = size > 0 ? new T[size] : nullptr;
    for (int i = 1; i < n_; ++i)
        v[i] = v[i - 1] + m_;
    for (int i = 0; i < n_; ++i)
        for (int j = 0; j < m_; ++j)
            v[i][j] = mat[i][j];
}

template <class T>
mymatrix<T> &mymatrix<T>::operator=(const mymatrix &mat)
{
    if (this != &mat)
    {
        if (n_ != mat.n_ || m_ != mat.m_)
        {
            if (v != nullptr)
            {
                delete[] v[0];
                delete[] v;
            }
            n_ = mat.n_;
            m_ = mat.m_;
            v = n_ > 0 ? new T *[n_] : nullptr;
            int size = n_ * m_;
            if (v)
                v[0] = size > 0 ? new T[size] : nullptr;
            for (int i = 1; i < n_; ++i)
                v[i] = v[i - 1] + m_;
        }
        for (int i = 0; i < n_; ++i)
            for (int j = 0; j < m_; ++j)
                v[i][j] = mat[i][j];
    }
    return *this;
}

template <class T>
T *mymatrix<T>::operator[](const int i)
{
    try
    {
        return v[i];
    }
    catch (const char *s)
    {
        throw std::logic_error("Matrix index is out of bounds.");
    }
}

template <class T>
const T *mymatrix<T>::operator[](const int i) const
{
    try
    {
        return v[i];
    }
    catch (const char *s)
    {
        throw std::logic_error("Matrix index is out of bounds.");
    }
}

template <class T>
mymatrix<T> mymatrix<T>::operator+(const mymatrix<T> &mat) const
{
    std::unique_ptr<mymatrix> temp(new mymatrix(*this));
    *temp += mat;
    return *temp;
}

template <class T>
mymatrix<T> mymatrix<T>::operator-(const mymatrix<T> &mat) const
{
    std::unique_ptr<mymatrix> temp(new mymatrix(*this));
    *temp -= mat;
    return *temp;
}

template <class T>
mymatrix<T> &mymatrix<T>::operator+=(const mymatrix<T> &mat)
{
    if (n_ != mat.n_ || m_ != mat.m_)
        throw std::logic_error("Sum possible only for the matrices of same size.");
    for (int i = 0; i < mat.n_; ++i)
        for (int j = 0; j < mat.m_; ++j)
            v[i][j] += mat[i][j];
    return *this;
}

template <class T>
mymatrix<T> &mymatrix<T>::operator-=(const mymatrix<T> &mat)
{
    if (n_ != mat.n_ || m_ != mat.m_)
        throw std::logic_error("Sum possible only for the matrices of same size.");
    for (int i = 0; i < mat.n_; ++i)
        for (int j = 0; j < mat.m_; ++j)
            v[i][j] -= mat[i][j];
    return *this;
}

template <typename T>
std::ostream &operator<<(std::ostream &os, const mymatrix<T> &mat)
{
    if (mat.v != nullptr)
    {
        os << "[[";
        for (int j = 0; j < mat.m_ - 1; ++j)
            os << mat[0][j] << ", ";
        os << mat[0][mat.m_ - 1] << "]\n";
        for (int i = 1; i < mat.n_ - 1; ++i)
        {
            os << " [";
            for (int j = 0; j < mat.m_ - 1; ++j)
                os << mat[i][j] << ", ";
            os << mat[i][mat.m_ - 1] << "]\n";
        }
        os << " [";
        for (int j = 0; j < mat.m_ - 1; ++j)
            os << mat[mat.n_ - 1][j] << ", ";
        os << mat[mat.n_ - 1][mat.m_ - 1] << "]";
        os << "]\n";
    }
    else
    {
        os << "Matrix undefined.\n";
    }
    return os;
}

template <class T>
int mymatrix<T>::nrows() const
{
    return n_;
}

template <class T>
int mymatrix<T>::ncols() const
{
    return m_;
}

template <class T>
void mymatrix<T>::resize(int newn, int newm)
{
    if (newn != n_ || newm != m_)
    {
        if (v != nullptr)
        {
            delete[] v[0];
            delete[] v;
        }
        n_ = newn;
        m_ = newm;
        v = n_ > 0 ? new T *[n_] : nullptr;
        int size = n_ * m_;
        if (v)
            v[0] = size > 0 ? new T[size] : nullptr;
        for (int i = 1; i < n_; ++i)
            v[i] = v[i - 1] + m_;
    }
}

template <class T>
void mymatrix<T>::assign(int newn, int newm, const T &a)
{
    if (newn != n_ || newm != m_)
    {
        if (v != nullptr)
        {
            delete[] v[0];
            delete[] v;
        }
        n_ = newn;
        m_ = newm;
        v = n_ > 0 ? new T *[n_] : nullptr;
        int size = n_ * m_;
        if (v)
            v[0] = size > 0 ? new T[size] : nullptr;
        for (int i = 1; i < n_; ++i)
            v[i] = v[i - 1] + m_;
    }
    for (int i = 0; i < n_; ++i)
        for (int j = 0; j < m_; ++j)
            v[i][j] = a;
}

template <class T>
void mymatrix<T>::identity()
{
    if (n_ != m_)
        throw std::runtime_error("Only square matrix can be identity.");
    for (int i = 0; i < n_; ++i)
    {
        for (int j = 0; j < m_; ++j)
            v[i][j] = 0;
        v[i][i] = 1;
    }
}

template <typename T>
mymatrix<T> mydot(const mymatrix<T> &mat1, const mymatrix<T> &mat2)
{
    mymatrix<T> result(mat1.n_, mat2.m_, 0);
    if (mat1.m_ != mat2.n_)
    {
        throw std::logic_error("Matrix product impossible.");
    }
    else
    {
        for (int i = 0; i < mat1.n_; ++i)
            for (int j = 0; j < mat2.m_; ++j)
                for (int k = 0; k < mat1.m_; ++k)
                    result[i][j] += mat1[i][k] * mat2[k][j];
    }
    return result;
}

template <typename T>
mymatrix<T> myouter(const mymatrix<T> &mat1, const mymatrix<T> &mat2)
{
    mymatrix<T> result(mat1.n_ * mat2.n_, mat1.m_ * mat2.m_, 0);
    for (int i = 0; i < result.n_; ++i)
        for (int j = 0; j < result.m_; ++j)
            result[i][j] += mat1[i / mat2.n_][j / mat2.m_] * mat2[i % mat2.n_][j % mat2.m_];

    return result;
}

template <typename T>
double mynorm(const mymatrix<T> &mat, const char *ord)
{
    if (mat.v == nullptr)
        throw std::logic_error("Norm of undefined matrix does not exist.");
    double result = 0;
    if (ord == "L2")
    {
        for (int i = 0; i < mat.n_; ++i)
            for (int j = 0; j < mat.m_; ++j)
                result += SQR(mat[i][j]);
    }
    else if (ord == "L1")
    {
        for (int i = 0; i < mat.n_; ++i)
            for (int j = 0; j < mat.m_; ++j)
                result += fabs(mat[i][j]);
    }
    else if (ord == "Linf")
    {
        for (int i = 0; i < mat.n_; ++i)
            for (int j = 0; j < mat.m_; ++j)
                if (fabs(mat[i][j]) > result)
                    result = fabs(mat[i][j]);
    }
    else
    {
        throw std::logic_error("Only L2, L1, Linf ord norm possible.");
    }
    return result;
}

template <class T>
mymatrix<T>::~mymatrix()
{
    if (v != nullptr)
    {
        delete[] v[0];
        delete[] v;
    }
}

typedef std::vector<std::vector<int>> initmatint;
typedef std::vector<std::vector<double>> initmatdoub;

typedef mymatrix<int> matint;
typedef mymatrix<double> matdoub;

#endif