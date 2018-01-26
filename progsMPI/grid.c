/*
This file is part of ftdock, a program for rigid-body protein-protein docking 
Copyright (C) 1997-2000 Gidon Moont

Biomolecular Modelling Laboratory
Imperial Cancer Research Fund
44 Lincoln's Inn Fields
London WC2A 3PX

+44 (0)20 7269 3348
http://www.bmm.icnet.uk/

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

*/

#include "structures.h"

void discretise_structure( struct Structure This_Structure , float grid_span , int grid_size , fftw_real *grid ) {

/************/

  /* Counters */

  int	residue , atom ;

  /* Co-ordinates */

  int	x , y , z ;
  int	steps , x_step , y_step , z_step ;

  float		x_centre , y_centre , z_centre ;

  /* Variables */

  float         distance , one_span ;

  /* create MPI atom & amino_acid type */

  MPI_Datatype MPI_Atom;//, MPI_Amino_Acid;

  #define atom_argc 6
  int atom_lengths[atom_argc] = {1,5,4,1,1,1};
  MPI_Datatype atom_types[atom_argc]
    = {MPI_INT, MPI_CHAR, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT};
  long atom_displs[atom_argc];
  atom_displs[0] = 0;
  atom_displs[1] = atom_displs[0] + sizeof(int);
  atom_displs[2] = atom_displs[1] + sizeof(char);
  atom_displs[3] = atom_displs[0] + sizeof(float);
  atom_displs[4] = atom_displs[0] + sizeof(float);
  atom_displs[5] = atom_displs[0] + sizeof(float);

  MPI_Type_create_struct( atom_argc, atom_lengths, atom_displs, atom_types, &MPI_Atom);
  MPI_Type_commit(&MPI_Atom);
  
  /*  #define amino_argc 7
  int amino_lengths[amino_argc] = {4,2,6,2,1,1,1};
  MPI_Datatype amino_types[amino_argc]
    = {MPI_CHAR, MPI_CHAR, MPI_CHAR, MPI_CHAR, MPI_INT, MPI_INT, MPI_Atom};
  long amino_displs[amino_argc];
  amino_displs[0] = 0;
  amino_displs[1] = amino_displs[0] + sizeof(char);
  amino_displs[2] = amino_displs[1] + sizeof(char);
  amino_displs[3] = amino_displs[2] + sizeof(char);
  amino_displs[4] = amino_displs[3] + sizeof(char);
  amino_displs[5] = amino_displs[4] + sizeof(int);
  amino_displs[6] = amino_displs[5] + sizeof(int);
  */


/************/

  one_span = grid_span / (float)grid_size ;

  distance = 1.8 ;

/************/
  
  /* create atom_array */
  int world_rank, world_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  int atom_num;
  struct Atom* atom_array;
  atom_num = 0;// could be computed in main
  for( residue = 1 ; residue <= This_Structure.length ; residue ++ ) {
    atom_num += This_Structure.Residue[residue].size;
  }
    
  atom_array = malloc(sizeof(struct Atom)*atom_num);

  int k = 0;
  for( residue = 1 ; residue <= This_Structure.length ; residue ++ ) {
    for( atom = 1 ; atom <= This_Structure.Residue[residue].size ; atom ++ ) {
      atom_array[k] = This_Structure.Residue[residue].Atom[atom];
      k++;
    }
  } 
  if (world_rank == 0) {
     for( x = 0 ; x < grid_size ; x ++ ) {
       for( y = 0 ; y < grid_size ; y ++ ) {
	 for( z = 0 ; z < grid_size ; z ++ ) {

	   grid[gaddress(x,y,z,grid_size)] = (fftw_real)0 ;

	 }
       }
     }
  }
/************/

  int workers = world_size - 1;
  
  //TODO: define my_atoms_lengths
  //TODO: define displs
  int atom_num_usual = atom_num / workers;
  int atom_num_remain= atom_num % workers;

  int *atom_num_each = (int*) malloc(sizeof(int) * world_size);
  int *displs_each   = (int*) malloc(sizeof(int) * world_size);
  char working       = workers;
  
  atom_num_each[0] = 0;
  displs_each[0]   = 0;
  for (int i=1; i < world_size; i++) {
    atom_num_each[i] = atom_num_usual + (i == 1 ? atom_num_remain : 0);
    displs_each[i]   = displs_each[i-1] + atom_num_each[i-1];
  }
  
  struct Atom *my_atoms = (struct Atom*) malloc(sizeof(struct Atom)*atom_num_each[world_rank]);

  MPI_Scatterv(atom_array, atom_num_each, displs_each, MPI_Atom,// sending args
  	       my_atoms, atom_num_each[world_rank], MPI_Atom, // receiving args
  	       0, MPI_COMM_WORLD); // other args

  steps = (int)( ( distance / one_span ) + 1.5 ) ;

  int sent;
  printf("%d atoms, %d working\n",atom_num_each[world_rank], working);
  if (world_rank == 0) {
      printf("%d\n", working);
    //root is the only one writing to grid
    MPI_Status status;
    while (working > 0) {
      MPI_Recv(&sent, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
      if (status.MPI_TAG == 0) {
	grid[sent] = (fftw_real)1 ;
      }
      if (status.MPI_TAG == 1) {
	working--;
      }
    }
  } else {
    //workers send gaddress 
    for (int i=0; i < atom_num_each[world_rank]; i++) {
      x = gord( my_atoms[i].coord[1] , grid_span , grid_size ) ;
      y = gord( my_atoms[i].coord[2] , grid_span , grid_size ) ;
      z = gord( my_atoms[i].coord[3] , grid_span , grid_size ) ;

      for( x_step = max( ( x - steps ) , 0 ) ;
	   x_step <= min( ( x + steps ) , ( grid_size - 1 ) ) ;
	   x_step ++ ) {

	x_centre  = gcentre( x_step , grid_span , grid_size ) ;

	for( y_step = max( ( y - steps ) , 0 ) ;
	     y_step <= min( ( y + steps ) , ( grid_size - 1 ) ) ;
	     y_step ++ ) {

	  y_centre  = gcentre( y_step , grid_span , grid_size ) ;

	  for( z_step = max( ( z - steps ) , 0 ) ;
	       z_step <= min( ( z + steps ) , ( grid_size - 1 ) )
		 ; z_step ++ ) {

	    z_centre  = gcentre( z_step , grid_span , grid_size ) ;

	    if( pythagoras( my_atoms[i].coord[1] ,
			    my_atoms[i].coord[2] ,
			    my_atoms[i].coord[3] ,
			    x_centre , y_centre , z_centre ) < distance ) {
	      sent = gaddress(x_step,y_step,z_step,grid_size) ;
	      MPI_Send(&sent, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
	    }
	  }
	}
    
      }
    }
  }
  MPI_Send(&sent, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);

/************/

  return ;

}



/************************/



void surface_grid( float grid_span , int grid_size , fftw_real *grid , float surface , float internal_value ) {


/************/

  /* Counters */

  int	x , y , z ;
  int	steps , x_step , y_step , z_step ;

  /* Variables */

  float		one_span ;

  int	at_surface ;

/************/

  one_span = grid_span / (float)grid_size ;

/************/

  /* Surface grid atoms */

  steps = (int)( ( surface / one_span ) + 1.5 ) ;

  for( x = 0 ; x < grid_size ; x ++ ) {
    for( y = 0 ; y < grid_size ; y ++ ) {
      for( z = 0 ; z < grid_size ; z ++ ) {

        if( (int)grid[gaddress(x,y,z,grid_size)] == 1 ) {

          at_surface = 0 ;

          for( x_step = max( x - steps , 0 ) ; x_step <= min( x + steps , grid_size - 1 ) ; x_step ++ ) {
            for( y_step = max( y - steps , 0 ) ; y_step <= min( y + steps , grid_size - 1 ) ; y_step ++ ) {
              for( z_step = max( z - steps , 0 ) ; z_step <= min( z + steps , grid_size - 1 ) ; z_step ++ ) {

                if( (int)grid[gaddress(x_step,y_step,z_step,grid_size)] == 0 ) {

                  if( ( (float)( ( ( x_step - x ) * ( x_step - x ) ) + ( ( y_step - y ) * ( y_step - y ) ) + ( ( z_step - z ) * ( z_step - z ) ) ) * one_span * one_span ) < ( surface * surface ) ) at_surface = 1 ;

                }

              }
            }
          }

          if( at_surface == 0 ) grid[gaddress(x,y,z,grid_size)] = (fftw_real)internal_value ;

        }

      }
    }
  }

/************/

  return ;

}
