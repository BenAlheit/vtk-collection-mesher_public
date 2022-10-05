//
// Created by alhei on 2022/09/09.
//

#ifndef VTK_COLLECTION_MESHER_SINGLEPARTICLE_H
#define VTK_COLLECTION_MESHER_SINGLEPARTICLE_H

#include "../Mesher.h"

template<unsigned int dim>
class SingleParticle{
public:
    SingleParticle();


private:

};

template<unsigned int dim>
SingleParticle<dim>::SingleParticle() {
    string fileName = "../examples/example_particle.vtk";
    string name = "in_prog";
    const unsigned int n_refines = 9;

    vtkNew <vtkPolyDataReader> reader;
    reader->SetFileName(fileName.c_str());
    reader->Update();
    vtkSmartPointer <vtkPolyData> polyData = reader->GetOutput();

    vector<vtkSmartPointer <vtkPolyData>> vtk_objects({polyData});
    vector<unsigned int> material_ids({1});

    Mesher<dim> mesher(name, vtk_objects, material_ids);
    mesher.write_mesh();
    for(unsigned int i=0; i<n_refines; i++){
        cout << "Refinement: " << i+1 << endl;
        mesher.refine();
        mesher.write_mesh();
    }

}

#endif //VTK_COLLECTION_MESHER_SINGLEPARTICLE_H
