#include <iostream>
#include "Mesher.h"
#include "PolycrystalIMPMesher.h"
#include "examples/SingleParticle.h"
#include "examples/Polycrystal.h"
#include "examples/PolycrystalIMPs.h"
using namespace std;

int main() {

//    SingleParticle<3> singleParticle = SingleParticle<3>();
//    Polycrystal<3> polycrystal = Polycrystal<3>();
    PolycrystalIMPs<3> polycrystal = PolycrystalIMPs<3>();
//    std::cout << "Hello, World!" << std::endl;
    return 0;
}
