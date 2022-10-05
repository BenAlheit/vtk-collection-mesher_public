//
// Created by alhei on 2022/09/09.
//

#ifndef VTK_COLLECTION_MESHER_POLYCRYSTAL_H
#define VTK_COLLECTION_MESHER_POLYCRYSTAL_H

//#include "../Mesher.h"
#include "../PolycrystalIMPMesher.h"
#include <vtkXMLPolyDataReader.h>

template<unsigned int dim>
class Polycrystal{
public:
    Polycrystal();


private:

};

template<unsigned int dim>
Polycrystal<dim>::Polycrystal() {
    string directory = "../examples/test-8g-0p/vtp/grain-";
    string name = "polycrystal";
    const unsigned int n_refines = 9;
    const unsigned int n_grains = 8;

    vector<vtkSmartPointer <vtkPolyData>> grains;

    for (unsigned int i = 0; i < n_grains; ++i) {
        vtkNew<vtkXMLPolyDataReader> reader;
        reader->SetFileName((directory+to_string(i+1)+".vtp").c_str());
        reader->Update();
        vtkSmartPointer <vtkPolyData> polyData = reader->GetOutput();
        grains.push_back(polyData);
    }

    PolycrystalIMPMesher<dim> mesher(name, grains);
    mesher.write_mesh();
    for(unsigned int i=0; i<n_refines; i++){
        cout << "Refinement: " << i+1 << endl;
        mesher.refine();
        mesher.write_mesh();
    }

}

#endif //VTK_COLLECTION_MESHER_POLYCRYSTAL_H
