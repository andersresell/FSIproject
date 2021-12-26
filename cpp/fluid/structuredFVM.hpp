//
// Created by anders on 12/25/21.
//

#ifndef FSIPROJECT_STRUCTUREDFVM_HPP
#define FSIPROJECT_STRUCTUREDFVM_HPP

#include <vector>
#include <iostream>

#define IS(i,j) (i*(ni+2) + j)
#define IV4(i,j) ()
static int ni,nj;
static int size = (ni+2)*(nj+2);

struct vec4{
    double u0, u1, u2, u3;
};

class FieldVec4{
    vec4* U;

public:
    FieldVec4();
    void operator*(double rhs);
    void operator*(FieldVec4& rhs);

    ~FieldVec4(){delete[] U;}
};




#endif //FSIPROJECT_STRUCTUREDFVM_HPP
