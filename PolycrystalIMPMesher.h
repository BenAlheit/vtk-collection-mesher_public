//
// Created by alhei on 2022/09/09.
//

#ifndef VTK_COLLECTION_MESHER_POLYCRYSTALIMPMESHER_H
#define VTK_COLLECTION_MESHER_POLYCRYSTALIMPMESHER_H

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
#include <vtkXMLPolyDataWriter.h>

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
#include <boost/archive/text_oarchive.hpp>

#include "Mesher.h"

using namespace std;
using namespace dealii;

bool almost_equals(const double &first, const double &second, const double &tol = 1e-10) {
    return fabs(first - second) < tol;
}

template<unsigned int dim>
class PolycrystalIMPMesher {
public:

    typedef pair <Point<dim>, Point<dim>> bounds;
    typedef pair <cell_identifier, cell_identifier> periodic_pair;

    PolycrystalIMPMesher(const string &name,
                         vector <vtkSmartPointer<vtkPolyData>> &grains,
                         const vector <vtkSmartPointer<vtkPolyData>> & particles = {},
                         const bool & grain_boundaries = true
    )
    : name (name)
    , grains(grains)
    , particles(particles)
    , grain_boundaries(grain_boundaries){
        system(("mkdir -p " + name).c_str());
        iota(range.begin(), range.end(), 0);
        double working_bounds[2 * dim], bounds[2 * dim];
        bool first = true;
        for (const auto &vtk_object: grains) {
            if (first) {
                vtk_object->GetCellsBounds(bounds);
            } else {
                vtk_object->GetCellsBounds(working_bounds);
                for (const auto &i: range) {
                    bounds[i * 2] = min(working_bounds[i * 2], bounds[i * 2]);
                    bounds[i * 2 + 1] = max(working_bounds[i * 2 + 1], bounds[i * 2 + 1]);
                }
            }

            first = false;
        }

        for (const auto &i: range) {
            p1[i] = bounds[i * 2];
            p2[i] = bounds[i * 2 + 1];
        }

        Point <dim> mid, diff;
        mid = (p1 + p2) / 2.;
        diff = (p2 - p1) / 2.;
        p2 = mid + 0.999 * diff;
        p1 = mid - 0.999 * diff;

        GridGenerator::hyper_rectangle(triangulation, p1, p2);

        unsigned int counter = 0;
        cell_identifier first_cell(0, 0);
        double p_bounds[6];
        pair < Point < dim > , Point < dim >> current_bounds = pair < Point < dim >, Point < dim
                >> (Point<dim>(), Point<dim>());
        for (const auto &vtk_object: particles) {
            for (vtkIdType i = 0; i < vtk_object->GetNumberOfPoints(); i++) {
                cell_to_vtk_points[first_cell].push_back(vtk_point_identifier(counter, i));
            }
            vtk_object->GetCellsBounds(p_bounds);
            for (const auto &i: range) {
                current_bounds.first[i] = p_bounds[2 * i];
                current_bounds.second[i] = p_bounds[2 * i + 1];
            }
            particle_bounds.push_back(current_bounds);
            counter++;
        }
        counter=0;
        for (const auto &vtk_object: particles) {
            string p_name = "particle-"+  to_string(counter) + ".vtp";
            string path = "./"+name+"/"+p_name;
            vtkNew<vtkXMLPolyDataWriter> writer;
            writer->SetFileName(path.c_str());
            writer->SetInputData(vtk_object);
            writer->Write();
            counter++;
        }



//        for (auto &face: triangulation.active_face_iterators()) {
//            for (const auto &i: range) {
//                if (almost_equals(face->center()[i], p1[i])) {
//                    face->set_boundary_id(i + 1);
//                    break;
//                }
//                if (almost_equals(face->center()[i], p2[i])) {
//                    face->set_boundary_id(dim + i + 1);
//                    break;
//                }
//            }
//        }

        update_vertex_to_grain_ids();
    };

    void refine();

    void write_mesh();

private:
    const string name;
    const bool grain_boundaries;
    vector <vtkSmartPointer<vtkPolyData>>
            grains;
    vector <vtkSmartPointer<vtkPolyData>>
            particles;

    vector <bounds> particle_bounds;

    vector <periodic_pair> periodic_pairs, previous_periodic_pairs;

    array<unsigned int, dim> range;

    unsigned int current_level = 0;
    const unsigned int particle_material_id = 1;

    unordered_map<unsigned int, unsigned int> vertex_to_grain_ids;

    unordered_map <cell_identifier, vector<vtk_point_identifier>, orderpairhash> cell_to_vtk_points;

    Triangulation <dim> triangulation;
    Point <dim> p1, p2;

    Point <dim> get_point(const vtk_point_identifier &id) {
        double x[3];
        Point <dim> out;
        particles.at(id.first)->GetPoints()->GetPoint(id.second, x);
        for (const auto &i: range) {
            out[i] = x[i];
        }
        return out;
    }

    bool in_bounds(Point <dim> pt, Point <dim> p1, Point <dim> p2) {
        bool out = true;
        for (const auto &i: range) {
            out = out && pt[i] <= p2[i];
            out = out && p1[i] <= pt[i];
        }
        return out;
    }

    bool in_bounds(Point <dim> pt, bounds _bounds) {
        bool out = true;
        for (const auto &i: range) {
            out = out && pt[i] <= _bounds.second[i];
            out = out && _bounds.first[i] <= pt[i];
        }
        return out;
    }

    void determine_materials();

    void update_vertex_to_grain_ids();

    void update_periodicity();

    void add_periodicity();

    bool is_periodic_pair(const CellAccessor <dim> &c1,
                          const CellAccessor <dim> &c2);

    void update_cell_to_vtk_points();

};

template<unsigned int dim>
void PolycrystalIMPMesher<dim>::refine() {

    CellAccessor <dim> current_cell;


    bool in_same_grain = true;
    for (const auto &cell: triangulation.active_cell_iterators_on_level(current_level)) {
        for (const auto i: cell->vertex_indices()) {
            for (const auto j: cell->vertex_indices()) {
//                in_same_grain = vertex_to_grain_ids.at(cell->vertex_iterator(i))
//                                == vertex_to_grain_ids.at(cell->vertex_iterator(j));
                in_same_grain = vertex_to_grain_ids.at(cell->vertex_index(i))
                                == vertex_to_grain_ids.at(cell->vertex_index(j));
                if (not in_same_grain)
                    break;
            }
            if (not in_same_grain)
                break;
        }
        if (not in_same_grain)
            cell->set_refine_flag();

        in_same_grain = true;
    }

    unsigned int level, cell_id;
    for (const auto &data: cell_to_vtk_points) {
        level = data.first.first;
        cell_id = data.first.second;
        if (level == current_level) {
            current_cell = CellAccessor<dim>(&triangulation, current_level, cell_id);
            if (not current_cell.refine_flag_set())
                current_cell.set_refine_flag();
        }
    }

    triangulation.execute_coarsening_and_refinement();
    current_level++;

    add_periodicity();

    update_cell_to_vtk_points();
    update_vertex_to_grain_ids();

    determine_materials();
    update_periodicity();
}

template<unsigned int dim>
void PolycrystalIMPMesher<dim>::write_mesh() {
    ofstream out("./" + name + "/" + name + "-" + to_string(current_level) + ".vtk");
    GridOut grid_out;
    grid_out.write_vtk(triangulation, out);

    ofstream output_file("./" + name + "/" + name + "-" + to_string(current_level) + ".dtri");
    {
        boost::archive::text_oarchive oa(output_file);
        triangulation.save(oa, 0);
    }
}

template<unsigned int dim>
void PolycrystalIMPMesher<dim>::determine_materials() {
    for (auto &cell: triangulation.active_cell_iterators_on_level(current_level)) {
        cell->set_material_id(0);
    }

    unordered_map<unsigned int, unsigned int> grain_to_count;
    unsigned int max_count, max_cell;
    max_count = 0;
    max_cell = 1000000; // TODO fix this horrific hack

    bool in_same_grain = true;
    if(grain_boundaries) {
        for (const auto &cell: triangulation.active_cell_iterators_on_level(current_level)) {
            for (const auto i: cell->vertex_indices()) {
                for (const auto j: cell->vertex_indices()) {
                    in_same_grain = vertex_to_grain_ids.at(cell->vertex_index(i))
                                    == vertex_to_grain_ids.at(cell->vertex_index(j));
                    if (not in_same_grain)
                        break;
                }
                if (not in_same_grain)
                    break;
            }
            if (in_same_grain)
                cell->set_material_id(vertex_to_grain_ids.at(cell->vertex_index(0)) + 2);

            in_same_grain = true;
        }
    }else {
//        for (const auto &cell: triangulation.active_cell_iterators_on_level(current_level)) {
////        for (const auto &cell: triangulation.active_cell_iterators()) {
//            grain_to_count.clear();
//            for (const auto i: cell->vertex_indices()) {
//                if (grain_to_count.find(vertex_to_grain_ids.at(cell->vertex_index(i))) == grain_to_count.end()) {
//                    grain_to_count[vertex_to_grain_ids.at(cell->vertex_index(i))] = 1;
//                } else {
//                    grain_to_count[vertex_to_grain_ids.at(cell->vertex_index(i))] += 1;
//                }
//            }
//
//            for (const auto &gc: grain_to_count) {
//                if (gc.second > max_count) {
//                    max_count = gc.second;
//                    max_cell = gc.first;
//                }
//            }
//
//            cell->set_material_id(vertex_to_grain_ids.at(max_cell) + 2);
//
//            max_count = 0;
//            max_cell = 1000000; // TODO fix this horrific hack
//
//        }
        for (const auto &cell: triangulation.active_cell_iterators_on_level(current_level)) {
            cell->set_material_id(vertex_to_grain_ids.at(cell->vertex_index(0)) + 2);
        }
    }

    unsigned int counter = 0;
    unsigned int cell_counter = 0;
    Point <dim> center, vtk_point;
    vector <cell_identifier> check_cells;
    cell_identifier cell_id;
    CellAccessor <dim> current_cell;

//    cell_id.first = current_level;
    for (const auto &particle: particles) {
        vtkSmartPointer <vtkPoints> check_centres = vtkSmartPointer<vtkPoints>::New();

//        for (const auto &cell: triangulation.active_cell_iterators_on_level(current_level)) {
        for (const auto &cell: triangulation.active_cell_iterators()) {
            center = cell->center();
            if (in_bounds(center, particle_bounds.at(counter))) {
                check_centres->InsertNextPoint(center[0],
                                               center[1],
                                               center[2]);
                cell_id.first = cell->level();
                cell_id.second = cell->index();
                check_cells.push_back(cell_id);
            }
        }

        vtkSmartPointer <vtkPolyData> pointsPolydata = vtkSmartPointer<vtkPolyData>::New();
        pointsPolydata->SetPoints(check_centres);

        vtkSmartPointer <vtkSelectEnclosedPoints> selectEnclosedPoints = vtkSmartPointer<vtkSelectEnclosedPoints>::New();
        selectEnclosedPoints->SetTolerance(1e-6);
        selectEnclosedPoints->SetSurfaceData(particle);
        selectEnclosedPoints->SetInputData(pointsPolydata);
        selectEnclosedPoints->Update();

        for (const auto &cell: check_cells) {
            if (selectEnclosedPoints->IsInside(cell_counter)) {
                current_cell = CellAccessor<dim>(&triangulation, cell.first, cell.second);
                current_cell.set_material_id(particle_material_id);
            }
            cell_counter++;
        }

        cell_counter = 0;

        check_cells.clear();
        counter++;
    }

}

template<unsigned int dim>
void PolycrystalIMPMesher<dim>::update_vertex_to_grain_ids() {
    vtkSmartPointer <vtkPoints> current_level_points = vtkSmartPointer<vtkPoints>::New();

    unordered_map<unsigned int, unsigned int> new_verts;

    unsigned int counter = 0;
    for (const auto &cell: triangulation.active_cell_iterators_on_level(current_level)) {
        for (const auto i: cell->vertex_indices()) {
            if (vertex_to_grain_ids.find(cell->vertex_index(i)) == vertex_to_grain_ids.end()) {
                current_level_points->InsertNextPoint(cell->vertex(i)[0],
                                                      cell->vertex(i)[1],
                                                      cell->vertex(i)[2]);
                new_verts[cell->vertex_index(i)] = counter;
                counter++;
            }
        }
    }

    vtkSmartPointer <vtkPolyData> pointsPolydata = vtkSmartPointer<vtkPolyData>::New();
    pointsPolydata->SetPoints(current_level_points);

    unsigned int obj_counter = 0;

    for (auto &vtk_object: grains) {
        vtkSmartPointer <vtkSelectEnclosedPoints> selectEnclosedPoints = vtkSmartPointer<vtkSelectEnclosedPoints>::New();
        selectEnclosedPoints->SetTolerance(1e-9);
        selectEnclosedPoints->SetSurfaceData(vtk_object);
        selectEnclosedPoints->SetInputData(pointsPolydata);
        selectEnclosedPoints->Update();

        for (const auto &vert: new_verts) {
            if (selectEnclosedPoints->IsInside(vert.second)) {
                vertex_to_grain_ids[vert.first] = obj_counter;
            }
        }
        obj_counter++;
    }
}

template<unsigned int dim>
void PolycrystalIMPMesher<dim>::update_periodicity() {
//    const CellAccessor<dim> c1, c2;
    cell_identifier c1_id, c2_id;
    c1_id.first = current_level;
    c2_id.first = current_level;

    typedef TriaAccessor <dim, dim, dim> cell_accessor;

    CellAccessor <dim> neg, pos, neg_child, pos_child;
    previous_periodic_pairs = periodic_pairs;
    periodic_pairs.clear();

//    TODO this method introduces duplicates which still works but isn't optimally efficient. Could improve by using a set and different typedef for periodic_pair where order doesn't matter
    if (previous_periodic_pairs.empty()) {
        for (const auto &c1: triangulation.active_cell_iterators_on_level(current_level)) {
            c1_id.second = c1->index();
            for (const auto &c2: triangulation.active_cell_iterators_on_level(current_level)) {
                if (c1->index() != c2->index())
                    if (is_periodic_pair(*c1, *c2)) {
                        c2_id.second = c2->index();

                        const periodic_pair pp(c1_id, c2_id);
                        periodic_pairs.push_back(pp);
                    }
            }
        }
    } else {
        for (const auto &parent_pp: previous_periodic_pairs) {
            if (parent_pp.first.first == current_level - 1) {
                neg = CellAccessor<dim>(&triangulation, parent_pp.first.first, parent_pp.first.second);
                pos = CellAccessor<dim>(&triangulation, parent_pp.second.first, parent_pp.second.second);
                for (unsigned int i = 0; i < neg.n_children(); ++i) {
                    neg_child = CellAccessor<dim>(&triangulation, current_level, neg.child_index(i));
                    for (unsigned int j = 0; j < pos.n_children(); ++j) {
                        pos_child = CellAccessor<dim>(&triangulation, current_level, pos.child_index(j));

                        if (is_periodic_pair(neg_child, pos_child)) {
                            c1_id.second = neg.child_index(i);
                            c2_id.second = pos.child_index(j);
                            const periodic_pair pp(c1_id, c2_id);
                            periodic_pairs.push_back(pp);
                        }
                    }
                }
            }
        }
    }
}

template<unsigned int dim>
void PolycrystalIMPMesher<dim>::add_periodicity() {
    unsigned int n_new_refines = 1;

    CellAccessor <dim> neg_cell, pos_cell;

    while (n_new_refines > 0) {
        n_new_refines = 0;
        for (const auto &pp: periodic_pairs) {
            if (pp.first.first == current_level - 1) {
                neg_cell = CellAccessor<dim>(&triangulation,
                                             current_level - 1,
                                             pp.first.second);
                pos_cell = CellAccessor<dim>(&triangulation,
                                             current_level - 1,
                                             pp.second.second);

                if (neg_cell.has_children() and not pos_cell.has_children()) {
                    pos_cell.set_refine_flag();
                    n_new_refines++;
                }
                if (not neg_cell.has_children() and pos_cell.has_children()) {
                    neg_cell.set_refine_flag();
                    n_new_refines++;
                }

            }
        }
        triangulation.execute_coarsening_and_refinement();
    }


}

template<unsigned int dim>
bool PolycrystalIMPMesher<dim>::is_periodic_pair(const CellAccessor <dim> &c1,
                                                 const CellAccessor <dim> &c2) {
    Tensor<1, dim> diff;
    unsigned int n_zeros = 0;

    for (const auto &f1: c1.face_iterators()) {
        if (f1->at_boundary()) {
            for (const auto &f2: c2.face_iterators()) {
                if (f2->at_boundary()) {
                    diff = f1->center() - f2->center();
                    for (const auto &i: range) {
                        n_zeros += (int) almost_equals(diff[i], 0);
                    }
                    if (n_zeros == 2)
                        return true;
                    n_zeros = 0;
                }
            }
        }
    }
//    double a = 0;
    return false;
}


template<unsigned int dim>
void PolycrystalIMPMesher<dim>::update_cell_to_vtk_points() {
    Point <dim> current_point;
    cell_identifier child_id, parent_id;

    child_id.first = current_level;
    parent_id.first = current_level - 1;
    for (const auto &cell: triangulation.active_cell_iterators_on_level(current_level)) {
        parent_id.second = cell->parent_index();

        if (cell_to_vtk_points.find(parent_id) != cell_to_vtk_points.end()) {
            const vector <vtk_point_identifier> &parent_points = cell_to_vtk_points.at(parent_id);

            child_id.second = cell->index();
            for (const auto &pt: parent_points) {
                current_point = get_point(pt);
                if (cell->point_inside(current_point)) {
                    cell_to_vtk_points[child_id].push_back(pt);
                }
            }
        }
    }
}

#endif //VTK_COLLECTION_MESHER_POLYCRYSTALIMPMESHER_H
