/*  ADMesh -- process triangulated solid meshes
 *  Copyright (C) 1995, 1996 Anthony D. Martin
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

#include <stdlib.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sstream>
#include "stl.h"

using namespace std;

#if !defined(SEEK_SET)
#define SEEK_SET 0
#define SEEK_CUR 1
#define SEEK_END 2
#endif

static void stl_put_little_int(FILE *fp, int value);
static void stl_put_little_float(FILE *fp, float value_in);
static void stl_initialize(stl *stl, char *file);
static void stl_allocate(stl *stl);
static void stl_read(stl *stl, int first_facet, int first);
static void stl_reallocate(stl *stl);
static int stl_get_little_int(FILE *fp);
static float stl_get_little_float(FILE *fp);

void stl::print_edges(FILE *file)
{
  int i;
  int edges_allocated;

  edges_allocated = stats.number_of_facets * 3;
  for(i = 0; i < edges_allocated; i++)
    {
      fprintf(file, "%d, %f, %f, %f, %f, %f, %f\n",
              edge_start[i].facet_number,
              edge_start[i].p1.x, edge_start[i].p1.y,
              edge_start[i].p1.z, edge_start[i].p2.x,
              edge_start[i].p2.y, edge_start[i].p2.z);
    }
}


void stl::stats_out(FILE *file, char *input_file)
{
  fprintf(file, "\n\
================= Results produced by ADMesh version 0.95 ================\n");
  fprintf(file, "\
Input file         : %s\n", input_file);
  if(stats.type == binary)
    {
      fprintf(file, "\
File type          : Binary STL file\n");
    }
  else
    {
      fprintf(file, "\
File type          : ASCII STL file\n");
    }      
  fprintf(file, "\
Header             : %s\n", stats.header);
  fprintf(file, "============== Size ==============\n");
  fprintf(file, "Min X = % f, Max X = % f\n", 
          stats.min.x, stats.max.x);
  fprintf(file, "Min Y = % f, Max Y = % f\n", 
          stats.min.y, stats.max.y);
  fprintf(file, "Min Z = % f, Max Z = % f\n", 
          stats.min.z, stats.max.z);
  
  fprintf(file, "\
========= Facet Status ========== Original ============ Final ====\n");
  fprintf(file, "\
Number of facets                 : %5d               %5d\n", 
          stats.original_num_facets, stats.number_of_facets);
  fprintf(file, "\
Facets with 1 disconnected edge  : %5d               %5d\n", 
          stats.facets_w_1_bad_edge, stats.connected_facets_2_edge -
          stats.connected_facets_3_edge);
  fprintf(file, "\
Facets with 2 disconnected edges : %5d               %5d\n", 
          stats.facets_w_2_bad_edge, stats.connected_facets_1_edge -
          stats.connected_facets_2_edge);
  fprintf(file, "\
Facets with 3 disconnected edges : %5d               %5d\n", 
          stats.facets_w_3_bad_edge, stats.number_of_facets -
          stats.connected_facets_1_edge);
  fprintf(file, "\
Total disconnected facets        : %5d               %5d\n", 
          stats.facets_w_1_bad_edge + stats.facets_w_2_bad_edge +
          stats.facets_w_3_bad_edge, stats.number_of_facets -
          stats.connected_facets_3_edge);
	  
  fprintf(file, 
"=== Processing Statistics ===     ===== Other Statistics =====\n");
  fprintf(file, "\
Number of parts       : %5d        Volume   : % f\n",
          stats.number_of_parts, stats.volume);
  fprintf(file, "\
Degenerate facets     : %5d\n", stats.degenerate_facets);
  fprintf(file, "\
Edges fixed           : %5d\n", stats.edges_fixed);
  fprintf(file, "\
Facets removed        : %5d\n", stats.facets_removed);
  fprintf(file, "\
Facets added          : %5d\n", stats.facets_added);
  fprintf(file, "\
Facets reversed       : %5d\n", stats.facets_reversed);
  fprintf(file, "\
Backwards edges       : %5d\n", stats.backwards_edges);
  fprintf(file, "\
Normals fixed         : %5d\n", stats.normals_fixed);
}

void stl::write_ascii(char *file, char *label)
{
  int       i;
  FILE      *fp;
  char      *error_msg;
  
  
  /* Open the file */
  fp = fopen(file, "w");
  if(fp == NULL)
    {
      error_msg = 
        (char *) malloc(81 + strlen(file)); /* Allow 80 chars+file size for message */
      sprintf(error_msg, "stl_write_ascii: Couldn't open %s for writing",
	      file);
      perror(error_msg);
      free(error_msg);
      exit(1);
    }
  
  fprintf(fp, "solid  %s\n", label);
  
  for(i = 0; i < stats.number_of_facets; i++)
    {
      fprintf(fp, "  facet normal % .8E % .8E % .8E\n",
              facet_start[i].normal.x, facet_start[i].normal.y,
              facet_start[i].normal.z);
      fprintf(fp, "    outer loop\n");
      fprintf(fp, "      vertex % .8E % .8E % .8E\n",
              facet_start[i].vertex[0].x, facet_start[i].vertex[0].y,
              facet_start[i].vertex[0].z);
      fprintf(fp, "      vertex % .8E % .8E % .8E\n",
              facet_start[i].vertex[1].x, facet_start[i].vertex[1].y,
              facet_start[i].vertex[1].z);
      fprintf(fp, "      vertex % .8E % .8E % .8E\n",
              facet_start[i].vertex[2].x, facet_start[i].vertex[2].y,
              facet_start[i].vertex[2].z);
      fprintf(fp, "    endloop\n");
      fprintf(fp, "  endfacet\n");
    }
  
  fprintf(fp, "endsolid  %s\n", label);
  
  fclose(fp);
}

void stl::print_neighbors(char *file)
{
  int i;
  FILE *fp;
  char *error_msg;

  /* Open the file */
  fp = fopen(file, "w");
  if(fp == NULL)
    {
      error_msg = 
        (char *) malloc(81 + strlen(file)); /* Allow 80 chars+file size for message */
      sprintf(error_msg, "stl_print_neighbors: Couldn't open %s for writing",
	      file);
      perror(error_msg);
      free(error_msg);
      exit(1);
    }

  for(i = 0; i < stats.number_of_facets; i++)
    {
      fprintf(fp, "%d, %d,%d, %d,%d, %d,%d\n",
	      i, 
              neighbors_start[i].neighbor[0],
              (int)neighbors_start[i].which_vertex_not[0],
              neighbors_start[i].neighbor[1],
              (int)neighbors_start[i].which_vertex_not[1],
              neighbors_start[i].neighbor[2],
              (int)neighbors_start[i].which_vertex_not[2]);
    }
}

static void stl_put_little_int(FILE *fp, int value_in)
{
  int new_value;
  union 
    {
      int  int_value;
      char char_value[4];
    } value;
  
  value.int_value = value_in;
  
  new_value  = value.char_value[0] & 0xFF;
  new_value |= (value.char_value[1] & 0xFF) << 0x08;
  new_value |= (value.char_value[2] & 0xFF) << 0x10;
  new_value |= (value.char_value[3] & 0xFF) << 0x18;
  fwrite(&new_value, sizeof(int), 1, fp);
}

static void stl_put_little_float(FILE *fp, float value_in)
{
  int new_value;
  union 
    {
      float float_value;
      char  char_value[4];
    } value;
  
  value.float_value = value_in;
  
  new_value  = value.char_value[0] & 0xFF;
  new_value |= (value.char_value[1] & 0xFF) << 0x08;
  new_value |= (value.char_value[2] & 0xFF) << 0x10;
  new_value |= (value.char_value[3] & 0xFF) << 0x18;
  fwrite(&new_value, sizeof(int), 1, fp);
}


void stl::write_binary(char *file, char *label)
{
  FILE      *fp;
  int       i;
  char      *error_msg;

  
  /* Open the file */
  fp = fopen(file, "w");
  if(fp == NULL)
    {
      error_msg = 
        (char *) malloc(81 + strlen(file)); /* Allow 80 chars+file size for message */
      sprintf(error_msg, "stl_write_binary: Couldn't open %s for writing",
	      file);
      perror(error_msg);
      free(error_msg);
      exit(1);
    }

  fprintf(fp, "%s", label);
  for(i = strlen(label); i < LABEL_SIZE; i++) putc(0, fp);

  fseek(fp, LABEL_SIZE, SEEK_SET);

  stl_put_little_int(fp, stats.number_of_facets);
  
  for(i = 0; i < stats.number_of_facets; i++)
    {
      stl_put_little_float(fp, facet_start[i].normal.x);
      stl_put_little_float(fp, facet_start[i].normal.y);
      stl_put_little_float(fp, facet_start[i].normal.z);
      stl_put_little_float(fp, facet_start[i].vertex[0].x);
      stl_put_little_float(fp, facet_start[i].vertex[0].y);
      stl_put_little_float(fp, facet_start[i].vertex[0].z);
      stl_put_little_float(fp, facet_start[i].vertex[1].x);
      stl_put_little_float(fp, facet_start[i].vertex[1].y);
      stl_put_little_float(fp, facet_start[i].vertex[1].z);
      stl_put_little_float(fp, facet_start[i].vertex[2].x);
      stl_put_little_float(fp, facet_start[i].vertex[2].y);
      stl_put_little_float(fp, facet_start[i].vertex[2].z);
      fputc(facet_start[i].extra[0], fp);
      fputc(facet_start[i].extra[1], fp);
    }
  
  fclose(fp);
}

void stl::write_vertex(int facet, int vertex)
{
  printf("  vertex %d/%d % .8E % .8E % .8E\n", vertex, facet,
         facet_start[facet].vertex[vertex].x,
         facet_start[facet].vertex[vertex].y,
         facet_start[facet].vertex[vertex].z);
}

void stl::write_facet(char *label, int facet)
{
  printf("facet (%d)/ %s\n", facet, label);
  write_vertex(facet, 0);
  write_vertex(facet, 1);
  write_vertex(facet, 2);
}

void stl::write_edge(char *label, stl_hash_edge edge)
{
  printf("edge (%d)/(%d) %s\n", edge.facet_number, edge.which_edge, label);
  if(edge.which_edge < 3)
    {
      write_vertex(edge.facet_number, edge.which_edge % 3);
      write_vertex(edge.facet_number, (edge.which_edge + 1) % 3);
    }
  else
    {
      write_vertex(edge.facet_number, (edge.which_edge + 1) % 3);
      write_vertex(edge.facet_number, edge.which_edge % 3);
    }
}

void stl::write_neighbor(int facet)
{
  printf("Neighbors %d: %d, %d, %d ;  %d, %d, %d\n", facet,
         neighbors_start[facet].neighbor[0],
         neighbors_start[facet].neighbor[1],
         neighbors_start[facet].neighbor[2],
         neighbors_start[facet].which_vertex_not[0],
         neighbors_start[facet].which_vertex_not[1],
         neighbors_start[facet].which_vertex_not[2]);
}

void stl::write_quad_object(char *file)
{
  FILE      *fp;
  int       i;
  int       j;
  char      *error_msg;
  stl_vertex connect_color;
  stl_vertex uncon_1_color;
  stl_vertex uncon_2_color;
  stl_vertex uncon_3_color;
  stl_vertex color;
  
  /* Open the file */
  fp = fopen(file, "w");
  if(fp == NULL)
    {
      error_msg = 
        (char *) malloc(81 + strlen(file)); /* Allow 80 chars+file size for message */
      sprintf(error_msg, "stl_write_quad_object: Couldn't open %s for writing",
	      file);
      perror(error_msg);
      free(error_msg);
      exit(1);
    }

  connect_color.x = 0.0;
  connect_color.y = 0.0;
  connect_color.z = 1.0;
  uncon_1_color.x = 0.0;
  uncon_1_color.y = 1.0;
  uncon_1_color.z = 0.0;
  uncon_2_color.x = 1.0;
  uncon_2_color.y = 1.0;
  uncon_2_color.z = 1.0;
  uncon_3_color.x = 1.0;
  uncon_3_color.y = 0.0;
  uncon_3_color.z = 0.0;

  fprintf(fp, "CQUAD\n");
  for(i = 0; i < stats.number_of_facets; i++)
    {
      j = ((neighbors_start[i].neighbor[0] == -1) +
           (neighbors_start[i].neighbor[1] == -1) +
           (neighbors_start[i].neighbor[2] == -1));
      if(j == 0)
	{
	  color = connect_color;
	}
      else if(j == 1)
	{
	  color = uncon_1_color;
	}
      else if(j == 2)
	{
	  color = uncon_2_color;
	}
      else
	{
	  color = uncon_3_color;
	}
      fprintf(fp, "%f %f %f    %1.1f %1.1f %1.1f 1\n", 
              facet_start[i].vertex[0].x,
              facet_start[i].vertex[0].y,
              facet_start[i].vertex[0].z, color.x, color.y, color.z);
      fprintf(fp, "%f %f %f    %1.1f %1.1f %1.1f 1\n",
              facet_start[i].vertex[1].x,
              facet_start[i].vertex[1].y,
              facet_start[i].vertex[1].z, color.x, color.y, color.z);
      fprintf(fp, "%f %f %f    %1.1f %1.1f %1.1f 1\n",
              facet_start[i].vertex[2].x,
              facet_start[i].vertex[2].y,
              facet_start[i].vertex[2].z, color.x, color.y, color.z);
      fprintf(fp, "%f %f %f    %1.1f %1.1f %1.1f 1\n",
              facet_start[i].vertex[2].x,
              facet_start[i].vertex[2].y,
              facet_start[i].vertex[2].z, color.x, color.y, color.z);
    }
  fclose(fp);
}
  
void stl::write_dxf(char *file, char *label)
{
  int       i;
  FILE      *fp;
  char      *error_msg;
  
  
  /* Open the file */
  fp = fopen(file, "w");
  if(fp == NULL)
    {
      error_msg = 
        (char *) malloc(81 + strlen(file)); /* Allow 80 chars+file size for message */
      sprintf(error_msg, "stl_write_ascii: Couldn't open %s for writing",
	      file);
      perror(error_msg);
      free(error_msg);
      exit(1);
    }
  
  fprintf(fp, "999\n%s\n", label);
  fprintf(fp, "0\nSECTION\n2\nHEADER\n0\nENDSEC\n");
  fprintf(fp, "0\nSECTION\n2\nTABLES\n0\nTABLE\n2\nLAYER\n70\n1\n\
0\nLAYER\n2\n0\n70\n0\n62\n7\n6\nCONTINUOUS\n0\nENDTAB\n0\nENDSEC\n");
  fprintf(fp, "0\nSECTION\n2\nBLOCKS\n0\nENDSEC\n");
  
  fprintf(fp, "0\nSECTION\n2\nENTITIES\n");

  for(i = 0; i < stats.number_of_facets; i++)
    {
      fprintf(fp, "0\n3DFACE\n8\n0\n");
      fprintf(fp, "10\n%f\n20\n%f\n30\n%f\n",
              facet_start[i].vertex[0].x, facet_start[i].vertex[0].y,
              facet_start[i].vertex[0].z);
      fprintf(fp, "11\n%f\n21\n%f\n31\n%f\n",
              facet_start[i].vertex[1].x, facet_start[i].vertex[1].y,
              facet_start[i].vertex[1].z);
      fprintf(fp, "12\n%f\n22\n%f\n32\n%f\n",
              facet_start[i].vertex[2].x, facet_start[i].vertex[2].y,
              facet_start[i].vertex[2].z);
      fprintf(fp, "13\n%f\n23\n%f\n33\n%f\n",
              facet_start[i].vertex[2].x, facet_start[i].vertex[2].y,
              facet_start[i].vertex[2].z);
    }
  
  fprintf(fp, "0\nENDSEC\n0\nEOF\n");
  
  fclose(fp);
}

void stl::open(char *file)
{
  stl_initialize(this, file);
  stl_allocate(this);
  stl_read(this, 0, 1);
  fclose(fp);
}

static int stl_get_little_int(FILE *fp)
{
  int value;
  value  =  fgetc(fp) & 0xFF;
  value |= (fgetc(fp) & 0xFF) << 0x08;
  value |= (fgetc(fp) & 0xFF) << 0x10;
  value |= (fgetc(fp) & 0xFF) << 0x18;
  return(value);
}

static float stl_get_little_float(FILE *fp)
{
  union
    {
      int   int_value;
      float float_value;
    } value;

  value.int_value  =  fgetc(fp) & 0xFF;
  value.int_value |= (fgetc(fp) & 0xFF) << 0x08;
  value.int_value |= (fgetc(fp) & 0xFF) << 0x10;
  value.int_value |= (fgetc(fp) & 0xFF) << 0x18;
  return(value.float_value);
}


static void stl_initialize(stl* stl, char *file)
{
  ulong           file_size;
  int            header_num_facets;
  int            num_facets;
  uint            i, j;
  unsigned char  chtest[128];
  int            num_lines = 1;
  char           *error_msg;

  stl->stats.degenerate_facets = 0;
  stl->stats.edges_fixed  = 0;
  stl->stats.facets_added = 0;
  stl->stats.facets_removed = 0;
  stl->stats.facets_reversed = 0;
  stl->stats.normals_fixed = 0;
  stl->stats.number_of_parts = 0;
  stl->stats.original_num_facets = 0;
  stl->stats.number_of_facets = 0;
  stl->stats.volume = -1.0;

  stl->neighbors_start = NULL;
  stl->facet_start = NULL;
  stl->v_indices = NULL;
  stl->v_shared = NULL;


  /* Open the file */
  stl->fp = fopen(file, "r");
  if(stl->fp == NULL)
    {
      error_msg =
        (char *) malloc(81 + strlen(file)); /* Allow 80 chars+file size for message */
      sprintf(error_msg, "stl_initialize: Couldn't open %s for reading",
              file);
      perror(error_msg);
      free(error_msg);
      exit(1);
    }
  /* Find size of file */
  fseek(stl->fp, 0, SEEK_END);
  file_size = ftell(stl->fp);

  /* Check for binary or ASCII file */
  fseek(stl->fp, HEADER_SIZE, SEEK_SET);
  fread(chtest, sizeof(chtest), 1, stl->fp);
  stl->stats.type = ascii;
  for(i = 0; i < sizeof(chtest); i++)
    {
      if(chtest[i] > 127)
        {
          stl->stats.type = binary;
          break;
        }
    }
  rewind(stl->fp);

  /* Get the header and the number of facets in the .STL file */
  /* If the .STL file is binary, then do the following */
  if(stl->stats.type == binary)
    {
      /* Test if the STL file has the right size  */
      if(((file_size - HEADER_SIZE) % SIZEOF_STL_FACET != 0)
         || (file_size < STL_MIN_FILE_SIZE))
        {
          fprintf(stderr, "The file %s has the wrong size.\n", file);
          exit(1);
        }
      num_facets = (file_size - HEADER_SIZE) / SIZEOF_STL_FACET;

      /* Read the header */
      fread(stl->stats.header, LABEL_SIZE, 1, stl->fp);
      stl->stats.header[80] = '\0';

      /* Read the int following the header.  This should contain # of facets */
      header_num_facets = stl_get_little_int(stl->fp);
      if(num_facets != header_num_facets)
        {
          fprintf(stderr,
          "Warning: File size doesn't match number of facets in the header\n");
        }
    }
  /* Otherwise, if the .STL file is ASCII, then do the following */
  else
    {
      /* Find the number of facets */
      j = 0;
      for(i = 0; i < file_size ; i++)
        {
          j++;
          if(getc(stl->fp) == '\n')
            {
              if(j > 4) /* don't count short lines */
                {
                  num_lines++;
                }
              j = 0;
            }
        }
      rewind(stl->fp);

      /* Get the header */
      for(i = 0;
          (i < 80) && (stl->stats.header[i] = getc(stl->fp)) != '\n'; i++);
      stl->stats.header[i] = '\0'; /* Lose the '\n' */
      stl->stats.header[80] = '\0';

      num_facets = num_lines / ASCII_LINES_PER_FACET;
    }
  stl->stats.number_of_facets += num_facets;
  stl->stats.original_num_facets = stl->stats.number_of_facets;
}

static void stl_allocate(stl* stl)
{
  /*  Allocate memory for the entire .STL file */
  stl->facet_start = (stl_facet*) calloc(stl->stats.number_of_facets,
                            sizeof(stl_facet));
  if(stl->facet_start == NULL) perror("stl_initialize");
  stl->stats.facets_malloced = stl->stats.number_of_facets;

  /* Allocate memory for the neighbors list */
  stl->neighbors_start =
    (stl_neighbors*) calloc(stl->stats.number_of_facets, sizeof(stl_neighbors));
  if(stl->facet_start == NULL) perror("stl_initialize");
}

void stl::open_merge(char *file)
{
  int first_facet;

  first_facet = stats.number_of_facets;
  stl_initialize(this, file);
  stl_reallocate(this);
  stl_read(this, first_facet, 0);
}

static void stl_reallocate(stl* stl)
{
  /*  Reallocate more memory for the .STL file(s) */
  stl->facet_start = (stl_facet*) realloc(stl->facet_start, stl->stats.number_of_facets *
                             sizeof(stl_facet));
  if(stl->facet_start == NULL) perror("stl_initialize");
  stl->stats.facets_malloced = stl->stats.number_of_facets;

  /* Reallocate more memory for the neighbors list */
  stl->neighbors_start =
    (stl_neighbors*)realloc(stl->neighbors_start, stl->stats.number_of_facets *
            sizeof(stl_neighbors));
  if(stl->facet_start == NULL) perror("stl_initialize");
}

static void stl_read(stl* stl, int first_facet, int first)
{
  stl_facet facet;
  int   i;
  float diff_x;
  float diff_y;
  float diff_z;
  float max_diff;


  if(stl->stats.type == binary)
    {
      fseek(stl->fp, HEADER_SIZE, SEEK_SET);
    }
  else
    {
      rewind(stl->fp);
      /* Skip the first line of the file */
      while(getc(stl->fp) != '\n');
    }

  for(i = first_facet; i < stl->stats.number_of_facets; i++)
    {
      if(stl->stats.type == binary)
        /* Read a single facet from a binary .STL file */
        {
          facet.normal.x = stl_get_little_float(stl->fp);
          facet.normal.y = stl_get_little_float(stl->fp);
          facet.normal.z = stl_get_little_float(stl->fp);
          facet.vertex[0].x = stl_get_little_float(stl->fp);
          facet.vertex[0].y = stl_get_little_float(stl->fp);
          facet.vertex[0].z = stl_get_little_float(stl->fp);
          facet.vertex[1].x = stl_get_little_float(stl->fp);
          facet.vertex[1].y = stl_get_little_float(stl->fp);
          facet.vertex[1].z = stl_get_little_float(stl->fp);
          facet.vertex[2].x = stl_get_little_float(stl->fp);
          facet.vertex[2].y = stl_get_little_float(stl->fp);
          facet.vertex[2].z = stl_get_little_float(stl->fp);
          facet.extra[0] = fgetc(stl->fp);
          facet.extra[1] = fgetc(stl->fp);
        }
      else
        /* Read a single facet from an ASCII .STL file */
        {
          fscanf(stl->fp, "%*s %*s %f %f %f\n", &facet.normal.x,
                 &facet.normal.y, &facet.normal.z);
          fscanf(stl->fp, "%*s %*s");
          fscanf(stl->fp, "%*s %f %f %f\n", &facet.vertex[0].x,
                 &facet.vertex[0].y,  &facet.vertex[0].z);
          fscanf(stl->fp, "%*s %f %f %f\n", &facet.vertex[1].x,
                 &facet.vertex[1].y,  &facet.vertex[1].z);
          fscanf(stl->fp, "%*s %f %f %f\n", &facet.vertex[2].x,
                 &facet.vertex[2].y,  &facet.vertex[2].z);
          fscanf(stl->fp, "%*s");
          fscanf(stl->fp, "%*s");
        }
      /* Write the facet into memory. */
      stl->facet_start[i] = facet;

      /* while we are going through all of the facets, let's find the  */
      /* maximum and minimum values for x, y, and z  */

      /* Initialize the max and min values the first time through*/
      if(first)
        {
          stl->stats.max.x = facet.vertex[0].x;
          stl->stats.min.x = facet.vertex[0].x;
          stl->stats.max.y = facet.vertex[0].y;
          stl->stats.min.y = facet.vertex[0].y;
          stl->stats.max.z = facet.vertex[0].z;
          stl->stats.min.z = facet.vertex[0].z;

          diff_x = ABS(facet.vertex[0].x - facet.vertex[1].x);
          diff_y = ABS(facet.vertex[0].y - facet.vertex[1].y);
          diff_z = ABS(facet.vertex[0].z - facet.vertex[1].z);
          max_diff = STL_MAX(diff_x, diff_y);
          max_diff = STL_MAX(diff_z, max_diff);
          stl->stats.shortest_edge = max_diff;

          first = 0;
        }
      /* now find the max and min values */
      stl->stats.max.x = STL_MAX(stl->stats.max.x, facet.vertex[0].x);
      stl->stats.min.x = STL_MIN(stl->stats.min.x, facet.vertex[0].x);
      stl->stats.max.y = STL_MAX(stl->stats.max.y, facet.vertex[0].y);
      stl->stats.min.y = STL_MIN(stl->stats.min.y, facet.vertex[0].y);
      stl->stats.max.z = STL_MAX(stl->stats.max.z, facet.vertex[0].z);
      stl->stats.min.z = STL_MIN(stl->stats.min.z, facet.vertex[0].z);

      stl->stats.max.x = STL_MAX(stl->stats.max.x, facet.vertex[1].x);
      stl->stats.min.x = STL_MIN(stl->stats.min.x, facet.vertex[1].x);
      stl->stats.max.y = STL_MAX(stl->stats.max.y, facet.vertex[1].y);
      stl->stats.min.y = STL_MIN(stl->stats.min.y, facet.vertex[1].y);
      stl->stats.max.z = STL_MAX(stl->stats.max.z, facet.vertex[1].z);
      stl->stats.min.z = STL_MIN(stl->stats.min.z, facet.vertex[1].z);

      stl->stats.max.x = STL_MAX(stl->stats.max.x, facet.vertex[2].x);
      stl->stats.min.x = STL_MIN(stl->stats.min.x, facet.vertex[2].x);
      stl->stats.max.y = STL_MAX(stl->stats.max.y, facet.vertex[2].y);
      stl->stats.min.y = STL_MIN(stl->stats.min.y, facet.vertex[2].y);
      stl->stats.max.z = STL_MAX(stl->stats.max.z, facet.vertex[2].z);
      stl->stats.min.z = STL_MIN(stl->stats.min.z, facet.vertex[2].z);
    }
  stl->stats.size.x = stl->stats.max.x - stl->stats.min.x;
  stl->stats.size.y = stl->stats.max.y - stl->stats.min.y;
  stl->stats.size.z = stl->stats.max.z - stl->stats.min.z;
  stl->stats.bounding_diameter =
    sqrt(stl->stats.size.x * stl->stats.size.x +
         stl->stats.size.y * stl->stats.size.y +
         stl->stats.size.z * stl->stats.size.z);
}


void stl::close()
{
    if(neighbors_start != NULL)
        free(neighbors_start);
    if(facet_start != NULL)
        free(facet_start);
    if(v_indices != NULL)
        free(v_indices);
    if(v_shared != NULL)
        free(v_shared);
}

void stl::generate_shared_vertices()
{
  int i;
  int j;
  int first_facet;
  int direction;
  int facet_num;
  int vnot;
  int next_edge;
  int pivot_vertex;
  int next_facet;
  int reversed;

  v_indices =
    (v_indices_struct*) calloc(stats.number_of_facets, sizeof(v_indices_struct));
  if(v_indices == NULL) perror("stl_generate_shared_vertices");
  v_shared =
    (stl_vertex*) calloc((stats.number_of_facets / 2), sizeof(stl_vertex));
  if(v_shared == NULL) perror("stl_generate_shared_vertices");
  stats.shared_malloced = stats.number_of_facets / 2;
  stats.shared_vertices = 0;

  for(i = 0; i < stats.number_of_facets; i++)
    {
      v_indices[i].vertex[0] = -1;
      v_indices[i].vertex[1] = -1;
      v_indices[i].vertex[2] = -1;
    }


  for(i = 0; i < stats.number_of_facets; i++)
    {
      first_facet = i;
      for(j = 0; j < 3; j++)
        {
          if(v_indices[i].vertex[j] != -1)
            {
              continue;
            }
          if(stats.shared_vertices == stats.shared_malloced)
            {
              stats.shared_malloced += 1024;
              v_shared = (stl_vertex*) realloc(v_shared,
                           stats.shared_malloced * sizeof(stl_vertex));
              if(v_shared == NULL) perror("stl_generate_shared_vertices");
            }

          v_shared[stats.shared_vertices] =
            facet_start[i].vertex[j];

          direction = 0;
          reversed = 0;
          facet_num = i;
          vnot = (j + 2) % 3;

          for(;;)
            {
              if(vnot > 2)
                {
                  if(direction == 0)
                    {
                      pivot_vertex = (vnot + 2) % 3;
                      next_edge = pivot_vertex;
                      direction = 1;
                    }
                  else
                    {
                      pivot_vertex = (vnot + 1) % 3;
                      next_edge = vnot % 3;
                      direction = 0;
                    }
                }
              else
                {
                  if(direction == 0)
                    {
                      pivot_vertex = (vnot + 1) % 3;
                      next_edge = vnot;
                    }
                  else
                    {
                      pivot_vertex = (vnot + 2) % 3;
                      next_edge = pivot_vertex;
                    }
                }
              v_indices[facet_num].vertex[pivot_vertex] =
                stats.shared_vertices;

              next_facet = neighbors_start[facet_num].neighbor[next_edge];
              if(next_facet == -1)
                {
                  if(reversed)
                    {
                      break;
                    }
                  else
                    {
                      direction = 1;
                      vnot = (j + 1) % 3;
                      reversed = 1;
                      facet_num = first_facet;
                    }
                }
              else if(next_facet != first_facet)
                {
                  vnot = neighbors_start[facet_num].
                    which_vertex_not[next_edge];
                  facet_num = next_facet;
                }
              else
                {
                  break;
                }
            }
          stats.shared_vertices += 1;
        }
    }
}

void stl::write_off(char *file)
{
  char      *error_msg;
  ofstream  outfile;

  outfile.open(file);

  if(outfile.fail())
    {
      error_msg =
        (char*) malloc(81 + strlen(file)); /* Allow 80 chars+file size for message */
      sprintf(error_msg, "write_off: Couldn't open %s for writing",
              file);
      perror(error_msg);
      free(error_msg);
      exit(1);
    }

  write_off(outfile);
  outfile.close();
}

void stl::write_off(ostream& stream)
{
    stream << "OFF\n";
    stream << stats.shared_vertices << " " << stats.number_of_facets << " 0" << endl;

    for(int i = 0; i < stats.shared_vertices; i++)
      {
        stream << "\t" << v_shared[i].x << " " << v_shared[i].y << " " << v_shared[i].z << endl;
      }
    for(int i = 0; i < stats.number_of_facets; i++)
      {
        fprintf(fp, "\t3 %d %d %d\n", v_indices[i].vertex[0],
                v_indices[i].vertex[1], v_indices[i].vertex[2]);
        stream << "\t3 " << v_indices[i].vertex[0] << " " << v_indices[i].vertex[1] << " " << v_indices[i].vertex[2] << endl;
      }
}

void stl::write_vrml(char *file)
{
  int i;
  FILE      *fp;
  char      *error_msg;


  /* Open the file */
  fp = fopen(file, "w");
  if(fp == NULL)
    {
      error_msg =
        (char*) malloc(81 + strlen(file)); /* Allow 80 chars+file size for message */
      sprintf(error_msg, "stl_write_ascii: Couldn't open %s for writing",
              file);
      perror(error_msg);
      free(error_msg);
      exit(1);
    }

  fprintf(fp, "#VRML V1.0 ascii\n\n");
  fprintf(fp, "Separator {\n");
  fprintf(fp, "\tDEF STLShape ShapeHints {\n");
  fprintf(fp, "\t\tvertexOrdering COUNTERCLOCKWISE\n");
  fprintf(fp, "\t\tfaceType CONVEX\n");
  fprintf(fp, "\t\tshapeType SOLID\n");
  fprintf(fp, "\t\tcreaseAngle 0.0\n");
  fprintf(fp, "\t}\n");
  fprintf(fp, "\tDEF STLModel Separator {\n");
  fprintf(fp, "\t\tDEF STLColor Material {\n");
  fprintf(fp, "\t\t\temissiveColor 0.700000 0.700000 0.000000\n");
  fprintf(fp, "\t\t}\n");
  fprintf(fp, "\t\tDEF STLVertices Coordinate3 {\n");
  fprintf(fp, "\t\t\tpoint [\n");

  for(i = 0; i < (stats.shared_vertices - 1); i++)
    {
      fprintf(fp, "\t\t\t\t%f %f %f,\n",
              v_shared[i].x, v_shared[i].y, v_shared[i].z);
    }
  fprintf(fp, "\t\t\t\t%f %f %f]\n",
          v_shared[i].x, v_shared[i].y, v_shared[i].z);
  fprintf(fp, "\t\t}\n");
  fprintf(fp, "\t\tDEF STLTriangles IndexedFaceSet {\n");
  fprintf(fp, "\t\t\tcoordIndex [\n");

  for(i = 0; i < (stats.number_of_facets - 1); i++)
    {
      fprintf(fp, "\t\t\t\t%d, %d, %d, -1,\n", v_indices[i].vertex[0],
              v_indices[i].vertex[1], v_indices[i].vertex[2]);
    }
  fprintf(fp, "\t\t\t\t%d, %d, %d, -1]\n", v_indices[i].vertex[0],
          v_indices[i].vertex[1], v_indices[i].vertex[2]);
  fprintf(fp, "\t\t}\n");
  fprintf(fp, "\t}\n");
  fprintf(fp, "}\n");
  fclose(fp);
}


Polyhedron stl::to_polyhedron()
{
    Polyhedron p;

    stringstream ss;

    write_off(ss);

    ss >> p;

    return p;
}
