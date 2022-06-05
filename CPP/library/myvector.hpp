#ifndef MYVECTOR_H
#define MYVECTOR_H

#include <istream>
#include <ostream>
#include <cmath>
#include <vector>

#include "mymatrix.hpp"

template <typename T>
class myvector;

template <typename T>
std::istream &operator>>(std::istream &is, myvector<T> &vec);

template <typename T>
std::ostream &operator<<(std::ostream &os, const myvector<T> &vec);

template <typename T>
double mydot(const myvector<T> &vec1, const myvector<T> &vec2);

template <typename T>
mymatrix<T> myouter(const myvector<T> &vec1, const myvector<T> &vec2);

template <typename T>
double mynorm(const myvector<T> &vec, const char *ord = "L2");

template <class T>
class myvector
{
private:
	int n_;
	T *v;

public:
	myvector();
	myvector(int n);
	myvector(int n, const T &a);
	myvector(int n, const std::vector<T> a);
	myvector(const myvector &vec);

	inline myvector &operator=(const myvector &vec);
	inline T &operator[](const int i);
	inline const T &operator[](const int i) const;
	inline myvector operator+(const myvector &vec) const;
	inline myvector operator-(const myvector &vec) const;
	inline myvector &operator+=(const myvector &vec);
	inline myvector &operator-=(const myvector &vec);
	friend std::istream &operator>><T>(std::istream &is, myvector<T> &vec);
	friend std::ostream &operator<<<T>(std::ostream &os, const myvector<T> &vec);

	inline int size() const;
	void resize(int newn);
	void assign(int newn, const T &a);
	friend double mydot<T>(const myvector<T> &vec1, const myvector<T> &vec2);
	friend mymatrix<T> myouter<T>(const myvector<T> &vec1, const myvector<T> &vec2);
	friend double mynorm<T>(const myvector<T> &vec, const char *ord);

	~myvector();
};

template <class T>
myvector<T>::myvector() : n_(0), v(nullptr) {}

template <class T>
myvector<T>::myvector(int n) : n_(n), v(n > 0 ? new T[n] : nullptr) {}

template <class T>
myvector<T>::myvector(int n, const T &a) : n_(n), v(n > 0 ? new T[n] : nullptr)
{
	for (int i = 0; i < n; ++i)
		v[i] = a;
}

template <class T>
myvector<T>::myvector(int n, const std::vector<T> a) : n_(n), v(n > 0 ? new T[n] : nullptr)
{
	for (int i = 0; i < n; ++i)
		v[i] = a[i];
}

template <class T>
myvector<T>::myvector(const myvector<T> &vec) : n_(vec.n_), v(n_ > 0 ? new T[n_] : nullptr)
{
	for (int i = 0; i < n_; ++i)
		v[i] = vec[i];
}

template <class T>
myvector<T> &myvector<T>::operator=(const myvector &vec)
{
	if (this != &vec)
	{
		if (n_ != vec.n_)
		{
			if (v != nullptr)
				delete[] v;
			n_ = vec.n_;
			v = n_ > 0 ? new T[n_] : nullptr;
		}
		for (int i = 0; i < n_; ++i)
			v[i] = vec[i];
	}
	return *this;
}

template <class T>
T &myvector<T>::operator[](const int i)
{
	try
	{
		return v[i];
	}
	catch (const char *s)
	{
		throw("Vector index is out of bounds.");
	}
}

template <class T>
const T &myvector<T>::operator[](const int i) const
{
	try
	{
		return v[i];
	}
	catch (const char *s)
	{
		throw("Vector index is out of bounds.");
	}
}

template <class T>
myvector<T> myvector<T>::operator+(const myvector &vec) const
{
	myvector temp(*this);
	temp += vec;
	return temp;
}

template <class T>
myvector<T> myvector<T>::operator-(const myvector &vec) const
{
	myvector temp(*this);
	temp += vec;
	return temp;
}

template <class T>
myvector<T> &myvector<T>::operator+=(const myvector &vec)
{
	if (n_ != vec.n_)
		throw("Sum possible only for the vectors of same size.");
	for (int i = 0; i < n_; ++i)
		v[i] += vec[i];
	return *this;
}

template <class T>
myvector<T> &myvector<T>::operator-=(const myvector &vec)
{
	if (n_ != vec.n_)
		throw("Sum possible only for the vectors of same size.");
	for (int i = 0; i < n_; ++i)
		v[i] -= vec[i];
	return *this;
}

template <typename T>
std::istream &operator>>(std::istream &is, myvector<T> &vec)
{
	if (vec.v != nullptr)
	{
		for (int i = 0; i < vec.n_; ++i)
			is >> vec[i];
	}
	else
	{
		throw("Vector is nullptr.");
	}
	return is;
}

template <typename T>
std::ostream &operator<<(std::ostream &os, const myvector<T> &vec)
{
	if (vec.v != nullptr)
	{
		os << "[";
		for (int i = 0; i < vec.n_ - 1; ++i)
			os << vec[i] << ", ";
		os << vec[vec.n_ - 1] << "]\n";
	}
	else
	{
		os << "Vector undefined.\n";
	}
	return os;
}

template <class T>
int myvector<T>::size() const
{
	return n_;
}

template <class T>
void myvector<T>::resize(int newn)
{
	if (newn != n_)
	{
		if (v != nullptr)
			delete[] v;
		n_ = newn;
		v = n_ > 0 ? new T[n_] : nullptr;
	}
}

template <class T>
void myvector<T>::assign(int newn, const T &a)
{
	if (newn != n_)
	{
		if (v != nullptr)
			delete[] v;
		n_ = newn;
		v = n_ > 0 ? new T[n_] : nullptr;
	}
	for (int i = 0; i < n_; ++i)
		v[i] = a;
}

template <typename T>
double mydot(const myvector<T> &vec1, const myvector<T> &vec2)
{
	double result = 0;
	if (vec1.n_ != vec2.n_)
		throw("Dot product possible only for vectors of the same size.");
	else
	{
		for (int i = 0; i < vec1.n_; ++i)
			result += vec1[i] * vec2[i];
	}
	return result;
}

template <typename T>
mymatrix<T> myouter(const myvector<T> &vec1, const myvector<T> &vec2)
{
	mymatrix<T> mat(vec1.n_, vec2.n_);
	for (int i = 0; i < vec1.n_; ++i)
		for (int j = 0; j < vec2.n_; ++j)
			mat[i][j] = vec1.v[i] * vec2.v[j];
	return mat;
}

template <typename T>
double mynorm(const myvector<T> &vec, const char *ord)
{
	if (vec.v == nullptr)
		throw("Norm of undefined vector does not exist.");
	double result = 0;
	if (ord == "L2")
	{
		for (int i = 0; i < vec.n_; ++i)
			result += vec[i] * vec[i];
		return sqrt(result);
	}
	else if (ord == "L1")
	{
		for (int i = 0; i < vec.n_; ++i)
			result += fabs(vec[i]);
		return result;
	}
	else if (ord == "Linf")
	{
		for (int i = 0; i < vec.n_; ++i)
			if (fabs(vec[i]) > result)
				result = fabs(vec[i]);
		return result;
	}
	else
	{
		throw("Only L2, L1, Linf ord norm possible.");
	}
}

template <class T>
myvector<T>::~myvector()
{
	if (v != nullptr)
		delete[] v;
}

typedef std::vector<int> initvecint;
typedef std::vector<double> initvecdoub;

typedef myvector<int> vecint;
typedef myvector<double> vecdoub;

#endif