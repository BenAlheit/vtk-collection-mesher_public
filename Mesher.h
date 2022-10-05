//
// Created by alhei on 2022/09/09.
//

#ifndef VTK_COLLECTION_MESHER_MESHER_H
#define VTK_COLLECTION_MESHER_MESHER_H

#include <vtkCylinderSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkProperty.h>
#include <vtkCamera.h>
#include <vtkSelectEnclosedPoints.h>
#include <vtkSmartPointer.h>
#include <vtkPolyDataReader.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkPolyDataNormals.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria_accessor.h>

#include <string>
#include <array>
#include <set>
#include <fstream>

#include <boost/functional/hash.hpp>

using namespace std;
using namespace dealii;

typedef pair<unsigned int, unsigned int> cell_identifier; // First -> level, Second -> cell index
typedef pair<unsigned int, unsigned int> vtk_point_identifier; // First -> vtk_object id, Second -> point id for vtk object

struct array_hash {
    template<class T1, long unsigned int n_elements>
    size_t operator()(const array<T1, n_elements> &p) const {
        return boost::hash_range(p.begin(), p.end());
    }
};

struct pairhash {
public:
    template<typename T, typename U>
    size_t operator()(const pair<T, U> &x) const {
        size_t seed = 0;
        boost::hash_combine(seed, x.first);
        boost::hash_combine(seed, x.second);
        return seed;
    }
};

struct orderpairhash {
public:
    template<typename T>
    size_t operator()(const pair<T, T> &x) const {
        array<T, 2> arr({x.first, x.second});
//        return array_hash<T, 2>(arr);
        return boost::hash_range(arr.begin(), arr.end());
    }
};

template<unsigned int dim>
class Mesher {
public:

    typedef pair<Point<dim>, Point<dim>> bounds;

    Mesher(const string &name,
           vector<vtkSmartPointer < vtkPolyData>>  & vtk_objects,
    const vector< unsigned int > & material_ids)
    : name (name)
    , vtk_objects(vtk_objects)
    , material_ids(material_ids)
    {
        system(("mkdir -p " + name).c_str());
        iota(range.begin(), range.end(), 0);
        double working_bounds[2 * dim], bounds[2 * dim];
        bool first = true;
        for (const auto &vtk_object: vtk_objects) {
            if (first) {
                vtk_object->GetCellsBounds(bounds);
            } else {
                vtk_object->GetCellsBounds(working_bounds);
                for (const auto &i: range) {
                    bounds[i] = min(working_bounds[i], bounds[i]);
                    bounds[i + dim] = max(working_bounds[i], bounds[i]);
                }
            }

            first = false;
        }

        for (const auto &i: range) {
            p1[i] = bounds[i*2];
            p2[i] = bounds[i*2 + 1];
        }

        GridGenerator::hyper_rectangle(triangulation, p1, p2);

        unsigned int counter = 0;
        cell_identifier first_cell(0, 0);
        for (const auto &vtk_object: vtk_objects) {
            for (vtkIdType i = 0; i < vtk_object->GetNumberOfPoints(); i++) {
                cell_to_vtk_points[first_cell].push_back(vtk_point_identifier(counter, i));
            }
            counter++;
        }

        for (auto &vtk_object: vtk_objects){

//            vtkSmartPointer <vtkPolyData> polydataCopy;
//            polydataCopy->DeepCopy(vtk_object);

            vtkNew <vtkPolyDataNormals> normalGenerator;
            normalGenerator->SetInputData(vtk_object);
            normalGenerator->ComputePointNormalsOn();
            normalGenerator->ComputeCellNormalsOff();
            normalGenerator->AutoOrientNormalsOn();
            normalGenerator->Update();
            vtk_object = normalGenerator->GetOutput();
            normals.push_back(vtk_object->GetPointData()->GetNormals());
        }
    };

    void refine();
    void write_mesh();

private:
    const string name;
    vector< vtkSmartPointer<vtkPolyData> > vtk_objects;
    const vector< unsigned int > material_ids;
    vector< vtkSmartPointer<vtkDataArray> > normals;
    vector< bounds > vtk_bounds;
    array<unsigned int, dim> range;

    unsigned int current_level = 0;

    unordered_map<cell_identifier, vector<vtk_point_identifier>, orderpairhash> cell_to_vtk_points;

    Triangulation<dim> triangulation;
    Point<dim> p1, p2;

    Point<dim> get_point(const vtk_point_identifier & id){
        double x[3];
        Point<dim> out;
        vtk_objects.at(id.first)->GetPoints()->GetPoint(id.second, x);
        for (const auto & i : range) {
            out[i] = x[i];
        }
        return out;
    }

    bool in_bounds(Point<dim> pt, Point<dim> p1, Point<dim> p2){
        bool out = true;
        for(const auto & i: range){
            out = out && pt[i] <= p2[i];
            out = out && p1[i] <= pt[i];
        }
        return out;
    }

    void determine_materials();
    void update_cell_to_vtk_points();
};

template<unsigned int dim>
void Mesher<dim>::refine() {

    CellAccessor<dim> current_cell;

    unsigned int level, cell_id;
    for (const auto &data: cell_to_vtk_points) {
        level = data.first.first;
        cell_id = data.first.second;
        if (level == current_level) {
            current_cell = CellAccessor<dim>(&triangulation, current_level, cell_id);
            if(not current_cell.refine_flag_set())
                current_cell.set_refine_flag();
        }
    }

    triangulation.execute_coarsening_and_refinement();
    current_level++;

    update_cell_to_vtk_points();
//    determine_materials();
}

template<unsigned int dim>
void Mesher<dim>::write_mesh() {
    ofstream out("./" + name + "/" + name + "-" + to_string(current_level) + ".vtk");
    GridOut grid_out;
    grid_out.write_vtk(triangulation, out);
}

template<unsigned int dim>
void Mesher<dim>::update_cell_to_vtk_points() {
    Point<dim> current_point;
    cell_identifier child_id, parent_id;

    child_id.first = current_level;
    parent_id.first = current_level - 1;
    for (const auto & cell: triangulation.active_cell_iterators_on_level(current_level)){
        parent_id.second = cell->parent_index();
        const vector<vtk_point_identifier> & parent_points = cell_to_vtk_points.at(parent_id);

        child_id.second = cell->index();
        for(const auto & pt: parent_points){
            current_point = get_point(pt);
            if(cell->point_inside(current_point)){
                cell_to_vtk_points[child_id].push_back(pt);
            }
        }
    }

}

template<unsigned int dim>
void Mesher<dim>::determine_materials() {
    Point<dim> center, vtk_point;
    Tensor<1, dim> diff;
    cell_identifier child_id, parent_id;

    child_id.first = current_level;
    parent_id.first = current_level - 1;
    bool contains_points;
    vector<vtk_point_identifier> close_points;
    double dot = 0;

    vector<cell_identifier> cells_without_points;
//    vector<Point<dim>> pointless_cell_centres;
    vtkSmartPointer <vtkPoints> pointless_cell_centres = vtkSmartPointer<vtkPoints>::New();

    for (const auto & cell: triangulation.active_cell_iterators_on_level(current_level)){
        center = cell->center();
        cell->set_material_id(0);
        child_id.second = cell->index();
        contains_points = cell_to_vtk_points.find(child_id) != cell_to_vtk_points.end();

        if (contains_points){
            close_points = cell_to_vtk_points.at(child_id);
        }else{
            parent_id.second = cell->parent_index();
            close_points = cell_to_vtk_points.at(parent_id);
        }

        if (contains_points) {
            for (const auto &pt: close_points) {
                vtk_point = get_point(pt);
                diff = center - vtk_point;
                for (const auto &i: range) {
                    dot += diff[i] * normals.at(pt.first)->GetComponent(pt.second, i);
                }
            }
        }else{
            cells_without_points.push_back(child_id);
//            pointless_cell_centres.push_back(center);
            pointless_cell_centres->InsertNextPoint((double) center[0], (double) center[1], (double) center[2]);
        }

        if(dot<0)
            cell->set_material_id(material_ids.at(close_points.at(0).first));

        dot = 0;
        contains_points = false;
    }

    CellAccessor<dim> current_cell;
    vtkSmartPointer <vtkPolyData> pointsPolydata = vtkSmartPointer<vtkPolyData>::New();
    pointsPolydata->SetPoints(pointless_cell_centres);

    unsigned int obj_counter = 0;
    unsigned int cell_counter = 0;

    for (auto &vtk_object: vtk_objects) {
        vtkSmartPointer <vtkSelectEnclosedPoints> selectEnclosedPoints = vtkSmartPointer<vtkSelectEnclosedPoints>::New();
        selectEnclosedPoints->SetSurfaceData(vtk_object);
        selectEnclosedPoints->SetInputData(pointsPolydata);
        selectEnclosedPoints->Update();

//        vector<cell_identifier> cells_without_points;
        for(const auto & cell: cells_without_points){
            if(selectEnclosedPoints->IsInside(cell_counter)){
                current_cell = CellAccessor<dim>(&triangulation, cell.first, cell.second);
                current_cell.set_material_id(material_ids.at(obj_counter));
            }
            cell_counter++;
        }

        cell_counter = 0;
        obj_counter++;
    }


}


#endif //VTK_COLLECTION_MESHER_MESHER_H
