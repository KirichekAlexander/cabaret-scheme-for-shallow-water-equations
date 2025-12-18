#include "aux_module.h"


//Реализация оператор для Vec3
InvC InvC::operator+(InvC const& other) const {
    return InvC{I + other.I,
                c};
}


InvC InvC::operator*(double k) const {
    return InvC{I * k,
                c};
}

InvC InvC::operator-(InvC const& other) const {
    return InvC{I - other.I,
                c};
}

InvC InvC::operator/(double k) const {
      return InvC{I / k,
                  c};
}


InvC InvC::operator*(InvC const& other) const {
    return InvC{I * other.I,
                c};
}


InvC operator*(double k, InvC const& vec) {
    return InvC{vec.I * k,
                vec.c};
}


Inv3 operator*(Vec3r const& vec1, Inv3 const& vec2) {
    return Inv3{InvC{vec1.a * vec2.a.I , vec2.a.c},
                InvC{vec1.b * vec2.b.I , vec2.b.c},
                InvC{vec1.c * vec2.c.I , vec2.c.c}};
}


Inv3 min_invariants(std::initializer_list<Inv3> list) {
    std::initializer_list<Inv3>::iterator it = list.begin();
    Inv3 res = *it++;
    for (; it != list.end(); ++it) {
        res.a = (res.a.I > (it->a).I ? it->a : res.a);
        res.b = (res.b.I > (it->b).I ? it->b : res.b);
        res.c = (res.c.I > (it->c).I ? it->c : res.c);
    }
    return res;
}


Inv3 max_invariants(std::initializer_list<Inv3> list) {
    std::initializer_list<Inv3>::iterator it = list.begin();
    Inv3 res = *it++;
    for (; it != list.end(); ++it) {
        res.a = (res.a.I < (it->a).I ? it->a : res.a);
        res.b = (res.b.I < (it->b).I ? it->b : res.b);
        res.c = (res.c.I < (it->c).I ? it->c : res.c);
    }
    return res;
}