/*
 * artio_mpi.h
 *
 *  Created on: Mar 6, 2012
 *      Author: Nick Gnedin
 */

#ifndef __ARTIO_MPI_H__
#define __ARTIO_MPI_H__

#include <mpi.h>

struct artio_context_struct {
	MPI_Comm comm;
};

#endif /* __ARTIO_MPI_H__ */
