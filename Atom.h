
#ifndef _ATOM_H_
#define _ATOM_H_

#include <string>
#include <sstream>
#include <iostream>

using namespace std;

class Atom {
private:
    double x, y, z;
public:
    Atom() = default;

    Atom(double x, double y, double z){
        this->x = x;
        this->y = y;
        this->z = z;
    }

    double getX() const {
        return x;
    }

    double getY() const {
        return y;
    }

    double getZ() const {
        return z;
    }

    // print particle in one string
    void printPositions(){
        cout<<"Atom positions: {" << x << "," << y << "," << z<<"}"<<endl;;
    }
};


#endif