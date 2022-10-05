# Vtk Collection Mesher

Applies a simple octree meshing algorithm (see chapter 8 of [this book](https://books.google.co.za/books/about/Automatic_Mesh_Generation.html?id=a6yPzQEACAAJ&redir_esc=y)) to [VTK](https://vtk.org/) poly-data objects using 
the [deal.ii](https://www.dealii.org/) finite element library.

# Functionality example: meshing a voronoi tesselation with complex geometries embedded in it

Consider the voronoi tesselation below that has complex geometries embedded in it.

![alt text](https://github.com/BenAlheit/vtk-collection-mesher_public/blob/master/imgs/voronois.png?raw=true)

This geometry is intended to be an approximation of the microstructure of a polycrystalline metal (the voronoi tesselation) that has intermetallic particles 
(the complex geometries, abbreviated as IMPs) embedded in it. The IMPs are created using image analysis of optical microscopy images (see this [this repo](https://github.com/BenAlheit/imp-image-analysis_public)
for details).

The meshing algorithm is incredibly simple:
1. Start with a coarse hex mesh (a single cell is fine);
2. Refine all the cells that intersect vtkPolyData objects;
3. Repeat step 2. until the mesh is a sufficiently close representation of the chosen geometry.

This mesh will, at a later stage, be used as a representative volume element (RVE) to conducted micromechanical modelling (see [this repo](https://github.com/BenAlheit/solver_public)) 
for which periodic boundary conditions are preferable. Hence, it is preferable to have a periodic mesh. This is ensured by refining elements that are "periodic pairs"; that is, 
if an element is part of a periodic pair and is refined then its partner must be refined too.

Ideally, there would be a final step that splits the final intersecting cells into a combination of hexahedral, tetrahedral, and wedge elements. However, 
the finite element library at hand only supports hexahedral elements. 

The resulting mesh for the geometry displayed above is presented below for increasing numbers of iterations.

![alt text](https://github.com/BenAlheit/vtk-collection-mesher_public/blob/master/imgs/meshes-together.png?raw=true)
