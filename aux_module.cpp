#include "aux_module.h"


//Реализация оператор для Vec3
Vec2 Vec2::operator+(Vec2 const& other) const {
    return Vec2{a + other.a,
                b};
}


Vec2 Vec2::operator*(double k) const {
    return Vec2{a * k,
                b};
}

Vec2 Vec2::operator-(Vec2 const& other) const {
    return Vec2{a - other.a,
                b};
}

Vec2 Vec2::operator/(double k) const {
      return Vec2{a / k,
                  b};
}


Vec2 Vec2::operator*(Vec2 const& other) const {
    return Vec2{a * other.a,
                b};
}


Vec2 operator*(double k, Vec2 const& vec) {
    return vec * k;
}


Vec3<Vec2> operator*(Vec3<double> const& vec1, Vec3<Vec2> const& vec2) {
    return Vec3<Vec2>{Vec2{vec1.a * vec2.a.a , vec2.a.b},
                      Vec2{vec1.b * vec2.b.a , vec2.b.b},
                      Vec2{vec1.c * vec2.c.a , vec2.c.b}};
}


Vec3<Vec2> min_invariants(std::initializer_list<Vec3<Vec2>> list) {
    std::initializer_list<Vec3<Vec2>>::iterator it = list.begin();
    Vec3<Vec2> res = *it++;
    for (; it != list.end(); ++it) {
        res.a = (res.a.a > (it->a).a ? it->a : res.a);
        res.b = (res.b.a > (it->b).a ? it->b : res.b);
        res.c = (res.c.a > (it->c).a ? it->c : res.c);
    }
    return res;
}


Vec3<Vec2> max_invariants(std::initializer_list<Vec3<Vec2>> list) {
    std::initializer_list<Vec3<Vec2>>::iterator it = list.begin();
    Vec3<Vec2> res = *it++;
    for (; it != list.end(); ++it) {
        res.a = (res.a.a < (it->a).a ? it->a : res.a);
        res.b = (res.b.a < (it->b).a ? it->b : res.b);
        res.c = (res.c.a < (it->c).a ? it->c : res.c);
    }
    return res;
}