#include <iostream>
#include <cmath>

// Dual number: a + e * b, e * e = 0, e != 0
template <typename T>
class Dual{
public:
    T real = T();
    T dual = T();
    
    Dual(T real_, T dual_) : real(real_), dual(dual_) {}
    Dual(T real_) : real(real_) {}

    auto operator+(const Dual<T>& y){
        return Dual<T>(this->real + y.real, this->dual + y.dual);
    }
    auto operator*(const Dual<T>& y){
        return Dual<T>(this->real * y.real, this->real * y.dual + this->dual * y.real);
    }
};

template <typename T>
Dual<T> forward(auto f, T primal, T tangent){
    Dual<T> input(primal, tangent);
    return f(input);
}

template <typename T>
T derivative(auto f, T real){
    auto result = forward(f, real, static_cast<T>(1));
    return result.dual;
}

// ---------- //

template <typename T>
T f(T x){
    return x * x;
}

template <typename T>
auto sin(Dual<T> x){
    return Dual<T>(std::sin(x.real), std::cos(x.real) * x.dual);
}

template<typename T>
auto sin(T x){
    return std::sin(x);
}

template <typename T>
T g(T x){
    return sin(x) + x;
}

int main(){
    auto x = 3.0;
    auto fx = f(x);
    auto dfx = derivative(f<Dual<double>>, x);
    std::cout << "f(" << x << ") = " << fx << "\ndf/dx(" << x << ") = " << dfx << "\n";

    x = 0.0;
    auto gx = g(x);
    auto dgx = derivative(g<Dual<double>>, x);
    std::cout << "g(" << x << ") = " << gx << "\ndg/dx(" << x << ") = " << dgx
              << "\n";
}