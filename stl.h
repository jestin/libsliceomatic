/*  ADMesh -- process triangulated solid meshes
 *  Copyright (C) 1995, 1996  Anthony D. Martin
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2, or (at your option)
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 *  
 *  Questions, comments, suggestions, etc to <amartin@engr.csulb.edu>
 */

#include <stdio.h>
#include <iostream>
#include <fstream>
#include "cgaldefs.h"

using namespace std;

#define STL_MAX(A,B) ((A)>(B)? (A):(B))
#define STL_MIN(A,B) ((A)<(B)? (A):(B))
#define ABS(X)  ((X) < 0 ? -(X) : (X))

#define LABEL_SIZE             80
#define NUM_FACET_SIZE         4
#define HEADER_SIZE            84
#define STL_MIN_FILE_SIZE      284
#define ASCII_LINES_PER_FACET  7
#define SIZEOF_EDGE_SORT       24

typedef struct 
{
  float x;
  float y;
  float z;
}stl_vertex;

typedef struct
{
  float x;
  float y;
  float z;
}stl_normal;

typedef char stl_extra[2];

typedef struct
{
  stl_normal normal;
  stl_vertex vertex[3];
  stl_extra  extra;
}stl_facet;
#define SIZEOF_STL_FACET       50

typedef enum {binary, ascii} stl_type;

typedef struct
{
  stl_vertex p1;
  stl_vertex p2;
  int        facet_number;
}stl_edge;

typedef struct stl_hash_edge
{
  unsigned       key[6];
  int            facet_number;
  int            which_edge;
  struct stl_hash_edge  *next;
}stl_hash_edge;

typedef struct
{
  int   neighbor[3];
  char  which_vertex_not[3];
}stl_neighbors;

typedef struct
{
  int   vertex[3];
}v_indices_struct;

typedef struct
{
  char          header[81];
  stl_type      type;
  int           number_of_facets;
  stl_vertex    max;
  stl_vertex    min;
  stl_vertex    size;
  float         bounding_diameter;
  float         shortest_edge;
  float         volume;
  unsigned      number_of_blocks;
  int           connected_edges;
  int           connected_facets_1_edge;
  int           connected_facets_2_edge;
  int           connected_facets_3_edge;
  int           facets_w_1_bad_edge;
  int           facets_w_2_bad_edge;
  int           facets_w_3_bad_edge;
  int           original_num_facets;
  int           edges_fixed;
  int           degenerate_facets;
  int           facets_removed;
  int           facets_added;
  int           facets_reversed;
  int           backwards_edges;
  int           normals_fixed;
  int           number_of_parts;
  int           malloced;
  int           freed;
  int           facets_malloced;
  int           collisions;
  int           shared_vertices;
  int           shared_malloced;
}stl_stats;

class stl
{
public:
    FILE          *fp;
    stl_facet     *facet_start;
    stl_edge      *edge_start;
    stl_hash_edge **heads;
    stl_hash_edge *tail;
    int           M;
    stl_neighbors *neighbors_start;
    v_indices_struct *v_indices;
    stl_vertex    *v_shared;
    stl_stats     stats;

    void open(char *file);
    void close();
    void stats_out(FILE *file, char *input_file);
    void print_edges(FILE *file);
    void print_neighbors(char *file);
    void write_ascii(char *file, char *label);
    void write_binary(char *file, char *label);
    void check_facets_exact();
    void check_facets_nearby(float tolerance);
    void remove_unconnected_facets();
    void write_vertex(int facet, int vertex);
    void write_facet(char *label, int facet);
    void write_edge(char *label, stl_hash_edge edge);
    void write_neighbor(int facet);
    void write_quad_object(char *file);
    void verify_neighbors();
    void fill_holes();
    void fix_normal_directions();
    void fix_normal_values();
    void reverse_all_facets();
    void translate(float x, float y, float z);
    void scale(float factor);
    void rotate_x(float angle);
    void rotate_y(float angle);
    void rotate_z(float angle);
    void mirror_xy();
    void mirror_yz();
    void mirror_xz();
    void open_merge(char *file);
    void generate_shared_vertices();
    void write_off(char *file);
    void write_off(ostream& stream);
    void write_dxf(char *file, char *label);
    void write_vrml(char *file);
    void calculate_normal(float normal[], stl_facet *facet);
    void normalize_vector(float v[]);
    void calculate_volume();

    Polyhedron to_polyhedron();
};
