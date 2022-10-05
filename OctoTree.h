/*
 * Ben Alheit, 12 Sept 2022:
 * This code was an idea that I was playing around with.
 * It isn't actually used in this project.
 *
 * */
#ifndef VTK_COLLECTION_MESHER_OCTOTREE_H
#define VTK_COLLECTION_MESHER_OCTOTREE_H



template<class Data>
class OctoNode {
public:
    OctoNode(const unsigned int &level,
             const unsigned int &level_id,
             const array<double, 3> &lower_bounds,
             const array<double, 3> &upper_bounds)
            : level(level), level_id(level_id), lower_bounds(lower_bounds), upper_bounds(upper_bounds) {};

protected:
    const unsigned int level;
    const unsigned int level_id;
    const array<double, 3> lower_bounds;
    const array<double, 3> upper_bounds;

private:
    array<OctoNode *, 8> next_nodes;

};

template<class Data>
class OctoLeaf : public OctoNode<Data> {
public:

    void make_leaves();

private:
    Data data;
};


template<class Data>
class OctoTree<Data> {
public:

    void refine();

private:
    OctoNode<Data> *head;

    unsigned int levels = 0;

};

class VtkCollection {
public:
    void get_point(const vtk_point_identifier & id, double pt[3]) const {
        vtk_objects.at(id.first)->GetPoints()->GetPoint(id.second, pt);
    };

private:
    const vector<vtkSmartPointer<vtkPolyData>> &vtk_objects;
};

class VtkCollectionSubset {
public:

    typedef vtk_point_identifier data_type;

    const vector<data_type> & get_iterator() const { return vtk_IDs;};

    bool in(const data_type &data,
            const array<double, 3> lower_bounds,
            const array<double, 3> upper_bounds) const {
        double x[3];
        parent_collection->get_point(data, x);
        return min(x[0], lower_bounds[0]) == lower_bounds[0] &&
               min(x[1], lower_bounds[1]) == lower_bounds[1] &&
               min(x[2], lower_bounds[2]) == lower_bounds[2] &&
               max(x[0], upper_bounds[0]) == upper_bounds[0] &&
               max(x[1], upper_bounds[1]) == upper_bounds[1] &&
               max(x[2], upper_bounds[2]) == upper_bounds[2];
    };

private:
    const VtkCollection *parent_collection;
    const vector<data_type> vtk_IDs;
};


#endif //VTK_COLLECTION_MESHER_OCTOTREE_H
