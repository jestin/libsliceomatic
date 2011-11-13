#ifndef CGALDEFS_H
#define CGALDEFS_H

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>


typedef CGAL::Simple_cartesian<double>  Kernel;
typedef Kernel::Point_3                 Point;
typedef Kernel::Vector_3                Vector;
typedef CGAL::Polyhedron_3<Kernel>      Polyhedron;
typedef Polyhedron::Vertex_iterator     Vertex_iterator;
typedef Polyhedron::Facet_iterator      Facet_iterator;
typedef Polyhedron::Halfedge_handle     Halfedge_handle;


#endif // CGALDEFS_H
