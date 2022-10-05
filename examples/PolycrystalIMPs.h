//
// Created by alhei on 2022/09/09.
//

#ifndef VTK_COLLECTION_MESHER_POLYCRYSTALIMPS_H
#define VTK_COLLECTION_MESHER_POLYCRYSTALIMPS_H

//#include "../Mesher.h"
#include "../PolycrystalIMPMesher.h"
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLPolyDataWriter.h>


template<unsigned int dim>
class PolycrystalIMPs {
public:
    PolycrystalIMPs();


private:
    array<unsigned int, dim> range;

    vtkSmartPointer <vtkPolyData> load_particle(string file,
                                                array<double, dim> length,
                                                array<double, dim> position);

};

template<unsigned int dim>
PolycrystalIMPs<dim>::PolycrystalIMPs() {
    iota(range.begin(), range.end(), 0);
    string directory = "../examples/test-8g-0p/vtp/grain-";
    string name = "polycrystal-imp-ngb";
//    string name = "polycrystal-ngb";
    const unsigned int n_refines = 7;
    const unsigned int n_grains = 8;

    vector<vtkSmartPointer < vtkPolyData>>
    grains;

    for (unsigned int i = 0; i < n_grains; ++i) {
        vtkNew <vtkXMLPolyDataReader> reader;
        reader->SetFileName((directory + to_string(i + 1) + ".vtp").c_str());
        reader->Update();
        vtkSmartPointer <vtkPolyData> polyData = reader->GetOutput();
        grains.push_back(polyData);
    }

    vector<vtkSmartPointer < vtkPolyData>>
    particles;

    string particle_name = "../examples/example_particle.vtk";
    array<double, dim> length({0.25, 0.25, 0.25});
    array<double, dim> position({0.25, 0.5, 0.5});

    vtkSmartPointer <vtkPolyData> polyData = load_particle(particle_name,
                                                           length,
                                                           position);

    array<double, dim> length2({0.25, 0.25, 0.4});
    array<double, dim> position2({0.6, 0.6, 0.35});
    vtkSmartPointer <vtkPolyData> polyData2 = load_particle(particle_name,
                                                            length2,
                                                            position2);

    particles.push_back(polyData);
    particles.push_back(polyData2);

    bool grain_boundaries = false;

    PolycrystalIMPMesher<dim> mesher(name, grains, particles, grain_boundaries);
    mesher.write_mesh();
    for (unsigned int i = 0; i < n_refines; i++) {
        cout << "Refinement: " << i + 1 << endl;
        mesher.refine();
        mesher.write_mesh();
    }

}

template<unsigned int dim>
vtkSmartPointer <vtkPolyData> PolycrystalIMPs<dim>::load_particle(string file,
                                                                  array<double, dim> length,
                                                                  array<double, dim> position) {
    vtkNew <vtkPolyDataReader> reader;
    reader->SetFileName(file.c_str());
    reader->Update();
    vtkSmartPointer <vtkPolyData> polyData = reader->GetOutput();

    double bounds[6];
    polyData->GetCellsBounds(bounds);
    array<double, dim> dist, center, scale, translate;
    for (const auto &i: range) {
        dist[i] = bounds[2 * i + 1] - bounds[2 * i];
        scale[i] = length[i] / dist[i];
    }

    vtkNew <vtkTransform> scaling;
    scaling->Scale(scale[0], scale[1], scale[2]);

    vtkNew <vtkTransformPolyDataFilter> scaleFilter;
    scaleFilter->SetInputData(polyData);
    scaleFilter->SetTransform(scaling);
    scaleFilter->Update();
    polyData = scaleFilter->GetOutput();

    polyData->GetCellsBounds(bounds);

    for (const auto &i: range) {
        center[i] = (bounds[2 * i + 1] + bounds[2 * i]) / 2.;
        translate[i] = position[i] - center[i];
    }

    vtkNew <vtkTransform> translation;
    translation->Translate(translate[0], translate[1], translate[2]);

    vtkNew <vtkTransformPolyDataFilter> transformFilter;
    transformFilter->SetInputData(polyData);
    transformFilter->SetTransform(translation);
    transformFilter->Update();
    polyData = transformFilter->GetOutput();

    polyData->GetCellsBounds(bounds);


    return polyData;
}

#endif //VTK_COLLECTION_MESHER_POLYCRYSTALIMPS_H
