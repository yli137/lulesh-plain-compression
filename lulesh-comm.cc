#include "lulesh.h"

#include <lz4.h>

#include <stdio.h>

// If no MPI, then this whole file is stubbed out
#if USE_MPI

#include <mpi.h>
#include <string.h>

/* Comm Routines */

#define ALLOW_UNPACKED_PLANE false
#define ALLOW_UNPACKED_ROW   false
#define ALLOW_UNPACKED_COL   false

int compress_lz4_buffer( const char *input_buffer, int input_size,
		char *output_buffer, int output_size )
{
	return LZ4_compress_default( input_buffer, output_buffer, input_size, output_size );
}


int decompress_lz4_buffer_default( const char *input_buffer, int input_size,
		char *output_buffer, int output_size )
{
	return LZ4_decompress_safe( input_buffer, output_buffer, input_size, output_size );
}

char *try_MPI_Isend( const void *buf, int count, MPI_Datatype type, int dest,
		int tag, MPI_Comm comm, MPI_Request *request )
{
	int size;
	MPI_Type_size( type, &size );
	size *= count;

	char *compress_buffer = (char*)malloc(size);
	int ret_size = compress_lz4_buffer( (const char*)buf, size, compress_buffer, size );

#if 0
	if( ret_size < size && size > 10000 ){
		MPI_Isend(compress_buffer, ret_size, MPI_BYTE, dest, tag, comm, request);
		return compress_buffer;
	}
#endif

	MPI_Isend(buf, count, type, dest, tag, comm, request);
	free(compress_buffer);
	return NULL;
}

void try_decompress( char *input_buffer, MPI_Status status )
{
	int count;
	MPI_Get_count(&status, MPI_CHAR, &count);

	char *output = (char*)malloc(count * 10);
	int output_size = count * 100;
	int ret = decompress_lz4_buffer_default( input_buffer, count, output, output_size );

	if( ret > count && ret > 0 ){
		printf("ret %7d count %7d\n", ret, count);
		memcpy(input_buffer, output, ret);
	}

	free(output);
}



// Initialize the dynamic request list
RequestList* init_request_list(int initial_capacity) {
    RequestList *list = (RequestList *)malloc(sizeof(RequestList));
    if (list == NULL) {
        fprintf(stderr, "Error: Could not allocate memory for request list\n");
        exit(EXIT_FAILURE);
    }
    list->entries = (RequestEntry *)malloc(initial_capacity * sizeof(RequestEntry));
    if (list->entries == NULL) {
        fprintf(stderr, "Error: Could not allocate memory for entries in request list\n");
        free(list); // Clean up previously allocated memory
        exit(EXIT_FAILURE);
    }
    list->size = 0;
    list->capacity = initial_capacity;
    return list;
}

// Function to append a new RequestEntry (MPI_Request and buffer address) to the list
void append_request(RequestList *list, MPI_Request *request, char *buffer) {
    // Check if the list needs to grow
    if (list->size >= list->capacity) {
        // Double the capacity
        list->capacity *= 2;
        list->entries = (RequestEntry *)realloc(list->entries, list->capacity * sizeof(RequestEntry));
        if (list->entries == NULL) {
            fprintf(stderr, "Error: Could not reallocate memory for request list\n");
            exit(EXIT_FAILURE);
        }
    }
    // Copy the MPI_Request into the list
    memcpy( list->entries[list->size].request, request, sizeof(MPI_Request) );
    list->entries[list->size].buffer = buffer;
    list->size++;
}

// Function to clean up the request list
void free_request_list(RequestList *list) {
    for (int i = 0; i < list->size; i++) {
        free(list->entries[i].buffer); // Free each buffer
    }
    free(list->entries);  // Free the array of RequestEntry
    free(list);           // Free the RequestList itself
}


/* doRecv flag only works with regular block structure */
RequestList *CommRecv(Domain& domain, Int_t msgType, Index_t xferFields,
		Index_t dx, Index_t dy, Index_t dz, bool doRecv, bool planeOnly) {

	if (domain.numRanks() == 1)
		return NULL;

	RequestList *request_list = init_request_list(26);

	/* post recieve buffers for all incoming messages */
	int myRank ;
	Index_t maxPlaneComm = xferFields * domain.maxPlaneSize() ;
	Index_t maxEdgeComm  = xferFields * domain.maxEdgeSize() ;
	Index_t pmsg = 0 ; /* plane comm msg */
	Index_t emsg = 0 ; /* edge comm msg */
	Index_t cmsg = 0 ; /* corner comm msg */
	MPI_Datatype baseType = ((sizeof(Real_t) == 4) ? MPI_FLOAT : MPI_DOUBLE) ;
	bool rowMin, rowMax, colMin, colMax, planeMin, planeMax ;

	/* assume communication to 6 neighbors by default */
	rowMin = rowMax = colMin = colMax = planeMin = planeMax = true ;

	if (domain.rowLoc() == 0) {
		rowMin = false ;
	}
	if (domain.rowLoc() == (domain.tp()-1)) {
		rowMax = false ;
	}
	if (domain.colLoc() == 0) {
		colMin = false ;
	}
	if (domain.colLoc() == (domain.tp()-1)) {
		colMax = false ;
	}
	if (domain.planeLoc() == 0) {
		planeMin = false ;
	}
	if (domain.planeLoc() == (domain.tp()-1)) {
		planeMax = false ;
	}

	for (Index_t i=0; i<26; ++i) {
		domain.recvRequest[i] = MPI_REQUEST_NULL ;
	}

	MPI_Comm_rank(MPI_COMM_WORLD, &myRank) ;

	/* post receives */

	/* receive data from neighboring domain faces */
	if (planeMin && doRecv) {
		/* contiguous memory */
		int fromRank = myRank - domain.tp()*domain.tp() ;
		int recvCount = dx * dy * xferFields ;
		MPI_Irecv(&domain.commDataRecv[pmsg * maxPlaneComm],
				recvCount, baseType, fromRank, msgType,
				MPI_COMM_WORLD, &domain.recvRequest[pmsg]);

		append_request(request_list, 
				&domain.recvRequest[pmsg], 
				(char*)&domain.commDataRecv[pmsg * maxPlaneComm]);

		++pmsg ;
	}
	if (planeMax) {
		/* contiguous memory */
		int fromRank = myRank + domain.tp()*domain.tp() ;
		int recvCount = dx * dy * xferFields ;
		MPI_Irecv(&domain.commDataRecv[pmsg * maxPlaneComm],
				recvCount, baseType, fromRank, msgType,
				MPI_COMM_WORLD, &domain.recvRequest[pmsg]) ;
		append_request(request_list, 
				&domain.recvRequest[pmsg], 
				(char*)&domain.commDataRecv[pmsg * maxPlaneComm]);
		++pmsg ;
	}
	if (rowMin && doRecv) {
		/* semi-contiguous memory */
		int fromRank = myRank - domain.tp() ;
		int recvCount = dx * dz * xferFields ;
		MPI_Irecv(&domain.commDataRecv[pmsg * maxPlaneComm],
				recvCount, baseType, fromRank, msgType,
				MPI_COMM_WORLD, &domain.recvRequest[pmsg]) ;
		append_request(request_list, 
				&domain.recvRequest[pmsg], 
				(char*)&domain.commDataRecv[pmsg * maxPlaneComm]);
		++pmsg ;
	}
	if (rowMax) {
		/* semi-contiguous memory */
		int fromRank = myRank + domain.tp() ;
		int recvCount = dx * dz * xferFields ;
		MPI_Irecv(&domain.commDataRecv[pmsg * maxPlaneComm],
				recvCount, baseType, fromRank, msgType,
				MPI_COMM_WORLD, &domain.recvRequest[pmsg]) ;
		append_request(request_list, 
				&domain.recvRequest[pmsg], 
				(char*)&domain.commDataRecv[pmsg * maxPlaneComm]);
		++pmsg ;
	}
	if (colMin && doRecv) {
		/* scattered memory */
		int fromRank = myRank - 1 ;
		int recvCount = dy * dz * xferFields ;
		MPI_Irecv(&domain.commDataRecv[pmsg * maxPlaneComm],
				recvCount, baseType, fromRank, msgType,
				MPI_COMM_WORLD, &domain.recvRequest[pmsg]) ;
		append_request(request_list, 
				&domain.recvRequest[pmsg], 
				(char*)&domain.commDataRecv[pmsg * maxPlaneComm]);
		++pmsg ;
	}
	if (colMax) {
		/* scattered memory */
		int fromRank = myRank + 1 ;
		int recvCount = dy * dz * xferFields ;
		MPI_Irecv(&domain.commDataRecv[pmsg * maxPlaneComm],
				recvCount, baseType, fromRank, msgType,
				MPI_COMM_WORLD, &domain.recvRequest[pmsg]) ;
		append_request(request_list, 
				&domain.recvRequest[pmsg], 
				(char*)&domain.commDataRecv[pmsg * maxPlaneComm]);
		++pmsg ;
	}

	if (!planeOnly) {
		/* receive data from domains connected only by an edge */
		if (rowMin && colMin && doRecv) {
			int fromRank = myRank - domain.tp() - 1 ;
			MPI_Irecv(&domain.commDataRecv[pmsg * maxPlaneComm +
					emsg * maxEdgeComm],
					dz * xferFields, baseType, fromRank, msgType,
					MPI_COMM_WORLD, &domain.recvRequest[pmsg+emsg]) ;
			append_request(request_list, 
					&domain.recvRequest[pmsg+emsg], 
					(char*)&domain.commDataRecv[pmsg * maxPlaneComm +
					emsg * maxEdgeComm]);
			++emsg ;
		}

		if (rowMin && planeMin && doRecv) {
			int fromRank = myRank - domain.tp()*domain.tp() - domain.tp() ;
			MPI_Irecv(&domain.commDataRecv[pmsg * maxPlaneComm +
					emsg * maxEdgeComm],
					dx * xferFields, baseType, fromRank, msgType,
					MPI_COMM_WORLD, &domain.recvRequest[pmsg+emsg]) ;
			append_request(request_list, 
					&domain.recvRequest[pmsg+emsg], 
					(char*)&domain.commDataRecv[pmsg * maxPlaneComm +
					emsg * maxEdgeComm]);
			++emsg ;
		}

		if (colMin && planeMin && doRecv) {
			int fromRank = myRank - domain.tp()*domain.tp() - 1 ;
			MPI_Irecv(&domain.commDataRecv[pmsg * maxPlaneComm +
					emsg * maxEdgeComm],
					dy * xferFields, baseType, fromRank, msgType,
					MPI_COMM_WORLD, &domain.recvRequest[pmsg+emsg]) ;
			append_request(request_list, 
					&domain.recvRequest[pmsg+emsg], 
					(char*)&domain.commDataRecv[pmsg * maxPlaneComm +
					emsg * maxEdgeComm]);
			++emsg ;
		}

		if (rowMax && colMax) {
			int fromRank = myRank + domain.tp() + 1 ;
			MPI_Irecv(&domain.commDataRecv[pmsg * maxPlaneComm +
					emsg * maxEdgeComm],
					dz * xferFields, baseType, fromRank, msgType,
					MPI_COMM_WORLD, &domain.recvRequest[pmsg+emsg]) ;
			append_request(request_list, 
					&domain.recvRequest[pmsg+emsg], 
					(char*)&domain.commDataRecv[pmsg * maxPlaneComm +
					emsg * maxEdgeComm]);
			++emsg ;
		}

		if (rowMax && planeMax) {
			int fromRank = myRank + domain.tp()*domain.tp() + domain.tp() ;
			MPI_Irecv(&domain.commDataRecv[pmsg * maxPlaneComm +
					emsg * maxEdgeComm],
					dx * xferFields, baseType, fromRank, msgType,
					MPI_COMM_WORLD, &domain.recvRequest[pmsg+emsg]) ;
			append_request(request_list, 
					&domain.recvRequest[pmsg+emsg], 
					(char*)&domain.commDataRecv[pmsg * maxPlaneComm +
					emsg * maxEdgeComm]);
			++emsg ;
		}

		if (colMax && planeMax) {
			int fromRank = myRank + domain.tp()*domain.tp() + 1 ;
			MPI_Irecv(&domain.commDataRecv[pmsg * maxPlaneComm +
					emsg * maxEdgeComm],
					dy * xferFields, baseType, fromRank, msgType,
					MPI_COMM_WORLD, &domain.recvRequest[pmsg+emsg]) ;
			append_request(request_list, 
					&domain.recvRequest[pmsg+emsg], 
					(char*)&domain.commDataRecv[pmsg * maxPlaneComm +
					emsg * maxEdgeComm]);
			++emsg ;
		}

		if (rowMax && colMin) {
			int fromRank = myRank + domain.tp() - 1 ;
			MPI_Irecv(&domain.commDataRecv[pmsg * maxPlaneComm +
					emsg * maxEdgeComm],
					dz * xferFields, baseType, fromRank, msgType,
					MPI_COMM_WORLD, &domain.recvRequest[pmsg+emsg]) ;
			append_request(request_list, 
					&domain.recvRequest[pmsg+emsg], 
					(char*)&domain.commDataRecv[pmsg * maxPlaneComm +
					emsg * maxEdgeComm]);
			++emsg ;
		}

		if (rowMin && planeMax) {
			int fromRank = myRank + domain.tp()*domain.tp() - domain.tp() ;
			MPI_Irecv(&domain.commDataRecv[pmsg * maxPlaneComm +
					emsg * maxEdgeComm],
					dx * xferFields, baseType, fromRank, msgType,
					MPI_COMM_WORLD, &domain.recvRequest[pmsg+emsg]) ;
			append_request(request_list, 
					&domain.recvRequest[pmsg+emsg], 
					(char*)&domain.commDataRecv[pmsg * maxPlaneComm +
					emsg * maxEdgeComm]);
			++emsg ;
		}

		if (colMin && planeMax) {
			int fromRank = myRank + domain.tp()*domain.tp() - 1 ;
			MPI_Irecv(&domain.commDataRecv[pmsg * maxPlaneComm +
					emsg * maxEdgeComm],
					dy * xferFields, baseType, fromRank, msgType,
					MPI_COMM_WORLD, &domain.recvRequest[pmsg+emsg]) ;
			append_request(request_list, 
					&domain.recvRequest[pmsg+emsg], 
					(char*)&domain.commDataRecv[pmsg * maxPlaneComm +
					emsg * maxEdgeComm]);
			++emsg ;
		}

		if (rowMin && colMax && doRecv) {
			int fromRank = myRank - domain.tp() + 1 ;
			MPI_Irecv(&domain.commDataRecv[pmsg * maxPlaneComm +
					emsg * maxEdgeComm],
					dz * xferFields, baseType, fromRank, msgType,
					MPI_COMM_WORLD, &domain.recvRequest[pmsg+emsg]) ;
			append_request(request_list, 
					&domain.recvRequest[pmsg+emsg], 
					(char*)&domain.commDataRecv[pmsg * maxPlaneComm +
					emsg * maxEdgeComm]);
			++emsg ;
		}

		if (rowMax && planeMin && doRecv) {
			int fromRank = myRank - domain.tp()*domain.tp() + domain.tp() ;
			MPI_Irecv(&domain.commDataRecv[pmsg * maxPlaneComm +
					emsg * maxEdgeComm],
					dx * xferFields, baseType, fromRank, msgType,
					MPI_COMM_WORLD, &domain.recvRequest[pmsg+emsg]) ;
			append_request(request_list, 
					&domain.recvRequest[pmsg+emsg], 
					(char*)&domain.commDataRecv[pmsg * maxPlaneComm +
					emsg * maxEdgeComm]);
			++emsg ;
		}

		if (colMax && planeMin && doRecv) {
			int fromRank = myRank - domain.tp()*domain.tp() + 1 ;
			MPI_Irecv(&domain.commDataRecv[pmsg * maxPlaneComm +
					emsg * maxEdgeComm],
					dy * xferFields, baseType, fromRank, msgType,
					MPI_COMM_WORLD, &domain.recvRequest[pmsg+emsg]) ;
			append_request(request_list, 
					&domain.recvRequest[pmsg+emsg], 
					(char*)&domain.commDataRecv[pmsg * maxPlaneComm +
					emsg * maxEdgeComm]);
			++emsg ;
		}

		/* receive data from domains connected only by a corner */
		if (rowMin && colMin && planeMin && doRecv) {
			/* corner at domain logical coord (0, 0, 0) */
			int fromRank = myRank - domain.tp()*domain.tp() - domain.tp() - 1 ;
			MPI_Irecv(&domain.commDataRecv[pmsg * maxPlaneComm +
					emsg * maxEdgeComm +
					cmsg * CACHE_COHERENCE_PAD_REAL],
					xferFields, baseType, fromRank, msgType,
					MPI_COMM_WORLD, &domain.recvRequest[pmsg+emsg+cmsg]) ;
			append_request(request_list, 
					&domain.recvRequest[pmsg+emsg+cmsg], 
					(char*)&domain.commDataRecv[pmsg * maxPlaneComm +
					emsg * maxEdgeComm +
					cmsg * CACHE_COHERENCE_PAD_REAL]);
			++cmsg ;
		}
		if (rowMin && colMin && planeMax) {
			/* corner at domain logical coord (0, 0, 1) */
			int fromRank = myRank + domain.tp()*domain.tp() - domain.tp() - 1 ;
			MPI_Irecv(&domain.commDataRecv[pmsg * maxPlaneComm +
					emsg * maxEdgeComm +
					cmsg * CACHE_COHERENCE_PAD_REAL],
					xferFields, baseType, fromRank, msgType,
					MPI_COMM_WORLD, &domain.recvRequest[pmsg+emsg+cmsg]) ;
			append_request(request_list, 
					&domain.recvRequest[pmsg+emsg+cmsg], 
					(char*)&domain.commDataRecv[pmsg * maxPlaneComm +
					emsg * maxEdgeComm +
					cmsg * CACHE_COHERENCE_PAD_REAL]);
			++cmsg ;
		}
		if (rowMin && colMax && planeMin && doRecv) {
			/* corner at domain logical coord (1, 0, 0) */
			int fromRank = myRank - domain.tp()*domain.tp() - domain.tp() + 1 ;
			MPI_Irecv(&domain.commDataRecv[pmsg * maxPlaneComm +
					emsg * maxEdgeComm +
					cmsg * CACHE_COHERENCE_PAD_REAL],
					xferFields, baseType, fromRank, msgType,
					MPI_COMM_WORLD, &domain.recvRequest[pmsg+emsg+cmsg]) ;
			append_request(request_list, 
					&domain.recvRequest[pmsg+emsg+cmsg], 
					(char*)&domain.commDataRecv[pmsg * maxPlaneComm +
					emsg * maxEdgeComm +
					cmsg * CACHE_COHERENCE_PAD_REAL]);
			++cmsg ;
		}
		if (rowMin && colMax && planeMax) {
			/* corner at domain logical coord (1, 0, 1) */
			int fromRank = myRank + domain.tp()*domain.tp() - domain.tp() + 1 ;
			MPI_Irecv(&domain.commDataRecv[pmsg * maxPlaneComm +
					emsg * maxEdgeComm +
					cmsg * CACHE_COHERENCE_PAD_REAL],
					xferFields, baseType, fromRank, msgType,
					MPI_COMM_WORLD, &domain.recvRequest[pmsg+emsg+cmsg]) ;
			append_request(request_list, 
					&domain.recvRequest[pmsg+emsg+cmsg], 
					(char*)&domain.commDataRecv[pmsg * maxPlaneComm +
					emsg * maxEdgeComm +
					cmsg * CACHE_COHERENCE_PAD_REAL]);
			++cmsg ;
		}
		if (rowMax && colMin && planeMin && doRecv) {
			/* corner at domain logical coord (0, 1, 0) */
			int fromRank = myRank - domain.tp()*domain.tp() + domain.tp() - 1 ;
			MPI_Irecv(&domain.commDataRecv[pmsg * maxPlaneComm +
					emsg * maxEdgeComm +
					cmsg * CACHE_COHERENCE_PAD_REAL],
					xferFields, baseType, fromRank, msgType,
					MPI_COMM_WORLD, &domain.recvRequest[pmsg+emsg+cmsg]) ;
			append_request(request_list, 
					&domain.recvRequest[pmsg+emsg+cmsg], 
					(char*)&domain.commDataRecv[pmsg * maxPlaneComm +
					emsg * maxEdgeComm +
					cmsg * CACHE_COHERENCE_PAD_REAL]);
			++cmsg ;
		}
		if (rowMax && colMin && planeMax) {
			/* corner at domain logical coord (0, 1, 1) */
			int fromRank = myRank + domain.tp()*domain.tp() + domain.tp() - 1 ;
			MPI_Irecv(&domain.commDataRecv[pmsg * maxPlaneComm +
					emsg * maxEdgeComm +
					cmsg * CACHE_COHERENCE_PAD_REAL],
					xferFields, baseType, fromRank, msgType,
					MPI_COMM_WORLD, &domain.recvRequest[pmsg+emsg+cmsg]) ;
			append_request(request_list, 
					&domain.recvRequest[pmsg+emsg+cmsg], 
					(char*)&domain.commDataRecv[pmsg * maxPlaneComm +
					emsg * maxEdgeComm +
					cmsg * CACHE_COHERENCE_PAD_REAL]);
			++cmsg ;
		}
		if (rowMax && colMax && planeMin && doRecv) {
			/* corner at domain logical coord (1, 1, 0) */
			int fromRank = myRank - domain.tp()*domain.tp() + domain.tp() + 1 ;
			MPI_Irecv(&domain.commDataRecv[pmsg * maxPlaneComm +
					emsg * maxEdgeComm +
					cmsg * CACHE_COHERENCE_PAD_REAL],
					xferFields, baseType, fromRank, msgType,
					MPI_COMM_WORLD, &domain.recvRequest[pmsg+emsg+cmsg]) ;
			append_request(request_list, 
					&domain.recvRequest[pmsg+emsg+cmsg], 
					(char*)&domain.commDataRecv[pmsg * maxPlaneComm +
					emsg * maxEdgeComm +
					cmsg * CACHE_COHERENCE_PAD_REAL]);
			++cmsg ;
		}
		if (rowMax && colMax && planeMax) {
			/* corner at domain logical coord (1, 1, 1) */
			int fromRank = myRank + domain.tp()*domain.tp() + domain.tp() + 1 ;
			MPI_Irecv(&domain.commDataRecv[pmsg * maxPlaneComm +
					emsg * maxEdgeComm +
					cmsg * CACHE_COHERENCE_PAD_REAL],
					xferFields, baseType, fromRank, msgType,
					MPI_COMM_WORLD, &domain.recvRequest[pmsg+emsg+cmsg]) ;
			append_request(request_list, 
					&domain.recvRequest[pmsg+emsg+cmsg], 
					(char*)&domain.commDataRecv[pmsg * maxPlaneComm +
					emsg * maxEdgeComm +
					cmsg * CACHE_COHERENCE_PAD_REAL]);
			++cmsg ;
		}
	}

	return request_list;
}

/******************************************/

void CommSend(Domain& domain, Int_t msgType,
		Index_t xferFields, Domain_member *fieldData,
		Index_t dx, Index_t dy, Index_t dz, bool doSend, bool planeOnly)
{
	char **addr = (char**)malloc(sizeof(char*) * 100);
	int ai = 0;
	for( int i = 0; i < 100; i++ ) addr[i] = NULL;


	if (domain.numRanks() == 1)
		return ;

	/* post recieve buffers for all incoming messages */
	int myRank ;
	Index_t maxPlaneComm = xferFields * domain.maxPlaneSize() ;
	Index_t maxEdgeComm  = xferFields * domain.maxEdgeSize() ;
	Index_t pmsg = 0 ; /* plane comm msg */
	Index_t emsg = 0 ; /* edge comm msg */
	Index_t cmsg = 0 ; /* corner comm msg */
	MPI_Datatype baseType = ((sizeof(Real_t) == 4) ? MPI_FLOAT : MPI_DOUBLE) ;
	MPI_Status status[26] ;
	Real_t *destAddr ;
	bool rowMin, rowMax, colMin, colMax, planeMin, planeMax ;
	/* assume communication to 6 neighbors by default */
	rowMin = rowMax = colMin = colMax = planeMin = planeMax = true ;
	if (domain.rowLoc() == 0) {
		rowMin = false ;
	}
	if (domain.rowLoc() == (domain.tp()-1)) {
		rowMax = false ;
	}
	if (domain.colLoc() == 0) {
		colMin = false ;
	}
	if (domain.colLoc() == (domain.tp()-1)) {
		colMax = false ;
	}
	if (domain.planeLoc() == 0) {
		planeMin = false ;
	}
	if (domain.planeLoc() == (domain.tp()-1)) {
		planeMax = false ;
	}

	for (Index_t i=0; i<26; ++i) {
		domain.sendRequest[i] = MPI_REQUEST_NULL ;
	}

	MPI_Comm_rank(MPI_COMM_WORLD, &myRank) ;

	/* post sends */

	if (planeMin | planeMax) {
		/* ASSUMING ONE DOMAIN PER RANK, CONSTANT BLOCK SIZE HERE */
		int sendCount = dx * dy ;

		if (planeMin) {
			destAddr = &domain.commDataSend[pmsg * maxPlaneComm] ;
			for (Index_t fi=0 ; fi<xferFields; ++fi) {
				Domain_member src = fieldData[fi] ;
				for (Index_t i=0; i<sendCount; ++i) {
					destAddr[i] = (domain.*src)(i) ;
				}
				destAddr += sendCount ;
			}
			destAddr -= xferFields*sendCount ;

			addr[ai++] = try_MPI_Isend(destAddr, xferFields*sendCount, baseType,
					myRank - domain.tp()*domain.tp(), msgType,
					MPI_COMM_WORLD, &domain.sendRequest[pmsg]) ;
			++pmsg ;
		}
		if (planeMax && doSend) {
			destAddr = &domain.commDataSend[pmsg * maxPlaneComm] ;
			for (Index_t fi=0 ; fi<xferFields; ++fi) {
				Domain_member src = fieldData[fi] ;
				for (Index_t i=0; i<sendCount; ++i) {
					destAddr[i] = (domain.*src)(dx*dy*(dz - 1) + i) ;
				}
				destAddr += sendCount ;
			}
			destAddr -= xferFields*sendCount ;

			addr[ai++] = try_MPI_Isend(destAddr, xferFields*sendCount, baseType,
					myRank + domain.tp()*domain.tp(), msgType,
					MPI_COMM_WORLD, &domain.sendRequest[pmsg]) ;
			++pmsg ;
		}
	}
	if (rowMin | rowMax) {
		/* ASSUMING ONE DOMAIN PER RANK, CONSTANT BLOCK SIZE HERE */
		int sendCount = dx * dz ;

		if (rowMin) {
			destAddr = &domain.commDataSend[pmsg * maxPlaneComm] ;
			for (Index_t fi=0; fi<xferFields; ++fi) {
				Domain_member src = fieldData[fi] ;
				for (Index_t i=0; i<dz; ++i) {
					for (Index_t j=0; j<dx; ++j) {
						destAddr[i*dx+j] = (domain.*src)(i*dx*dy + j) ;
					}
				}
				destAddr += sendCount ;
			}
			destAddr -= xferFields*sendCount ;

			addr[ai++] = try_MPI_Isend(destAddr, xferFields*sendCount, baseType,
					myRank - domain.tp(), msgType,
					MPI_COMM_WORLD, &domain.sendRequest[pmsg]) ;
			++pmsg ;
		}
		if (rowMax && doSend) {
			destAddr = &domain.commDataSend[pmsg * maxPlaneComm] ;
			for (Index_t fi=0; fi<xferFields; ++fi) {
				Domain_member src = fieldData[fi] ;
				for (Index_t i=0; i<dz; ++i) {
					for (Index_t j=0; j<dx; ++j) {
						destAddr[i*dx+j] = (domain.*src)(dx*(dy - 1) + i*dx*dy + j) ;
					}
				}
				destAddr += sendCount ;
			}
			destAddr -= xferFields*sendCount ;

			addr[ai++] = try_MPI_Isend(destAddr, xferFields*sendCount, baseType,
					myRank + domain.tp(), msgType,
					MPI_COMM_WORLD, &domain.sendRequest[pmsg]) ;
			++pmsg ;
		}
	}
	if (colMin | colMax) {
		/* ASSUMING ONE DOMAIN PER RANK, CONSTANT BLOCK SIZE HERE */
		int sendCount = dy * dz ;

		if (colMin) {
			destAddr = &domain.commDataSend[pmsg * maxPlaneComm] ;
			for (Index_t fi=0; fi<xferFields; ++fi) {
				Domain_member src = fieldData[fi] ;
				for (Index_t i=0; i<dz; ++i) {
					for (Index_t j=0; j<dy; ++j) {
						destAddr[i*dy + j] = (domain.*src)(i*dx*dy + j*dx) ;
					}
				}
				destAddr += sendCount ;
			}
			destAddr -= xferFields*sendCount ;

			addr[ai++] = try_MPI_Isend(destAddr, xferFields*sendCount, baseType,
					myRank - 1, msgType,
					MPI_COMM_WORLD, &domain.sendRequest[pmsg]) ;
			++pmsg ;
		}
		if (colMax && doSend) {
			destAddr = &domain.commDataSend[pmsg * maxPlaneComm] ;
			for (Index_t fi=0; fi<xferFields; ++fi) {
				Domain_member src = fieldData[fi] ;
				for (Index_t i=0; i<dz; ++i) {
					for (Index_t j=0; j<dy; ++j) {
						destAddr[i*dy + j] = (domain.*src)(dx - 1 + i*dx*dy + j*dx) ;
					}
				}
				destAddr += sendCount ;
			}
			destAddr -= xferFields*sendCount ;

			addr[ai++] = try_MPI_Isend(destAddr, xferFields*sendCount, baseType,
					myRank + 1, msgType,
					MPI_COMM_WORLD, &domain.sendRequest[pmsg]) ;
			++pmsg ;
		}
	}

	if (!planeOnly) {
		if (rowMin && colMin) {
			int toRank = myRank - domain.tp() - 1 ;
			destAddr = &domain.commDataSend[pmsg * maxPlaneComm +
				emsg * maxEdgeComm] ;
			for (Index_t fi=0; fi<xferFields; ++fi) {
				Domain_member src = fieldData[fi] ;
				for (Index_t i=0; i<dz; ++i) {
					destAddr[i] = (domain.*src)(i*dx*dy) ;
				}
				destAddr += dz ;
			}
			destAddr -= xferFields*dz ;
			addr[ai++] = try_MPI_Isend(destAddr, xferFields*dz, baseType, toRank, msgType,
					MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg]) ;
			++emsg ;
		}

		if (rowMin && planeMin) {
			int toRank = myRank - domain.tp()*domain.tp() - domain.tp() ;
			destAddr = &domain.commDataSend[pmsg * maxPlaneComm +
				emsg * maxEdgeComm] ;
			for (Index_t fi=0; fi<xferFields; ++fi) {
				Domain_member src = fieldData[fi] ;
				for (Index_t i=0; i<dx; ++i) {
					destAddr[i] = (domain.*src)(i) ;
				}
				destAddr += dx ;
			}
			destAddr -= xferFields*dx ;
			addr[ai++] = try_MPI_Isend(destAddr, xferFields*dx, baseType, toRank, msgType,
					MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg]) ;
			++emsg ;
		}

		if (colMin && planeMin) {
			int toRank = myRank - domain.tp()*domain.tp() - 1 ;
			destAddr = &domain.commDataSend[pmsg * maxPlaneComm +
				emsg * maxEdgeComm] ;
			for (Index_t fi=0; fi<xferFields; ++fi) {
				Domain_member src = fieldData[fi] ;
				for (Index_t i=0; i<dy; ++i) {
					destAddr[i] = (domain.*src)(i*dx) ;
				}
				destAddr += dy ;
			}
			destAddr -= xferFields*dy ;
			addr[ai++] = try_MPI_Isend(destAddr, xferFields*dy, baseType, toRank, msgType,
					MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg]) ;
			++emsg ;
		}

		if (rowMax && colMax && doSend) {
			int toRank = myRank + domain.tp() + 1 ;
			destAddr = &domain.commDataSend[pmsg * maxPlaneComm +
				emsg * maxEdgeComm] ;
			for (Index_t fi=0; fi<xferFields; ++fi) {
				Domain_member src = fieldData[fi] ;
				for (Index_t i=0; i<dz; ++i) {
					destAddr[i] = (domain.*src)(dx*dy - 1 + i*dx*dy) ;
				}
				destAddr += dz ;
			}
			destAddr -= xferFields*dz ;
			addr[ai++] = try_MPI_Isend(destAddr, xferFields*dz, baseType, toRank, msgType,
					MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg]) ;
			++emsg ;
		}

		if (rowMax && planeMax && doSend) {
			int toRank = myRank + domain.tp()*domain.tp() + domain.tp() ;
			destAddr = &domain.commDataSend[pmsg * maxPlaneComm +
				emsg * maxEdgeComm] ;
			for (Index_t fi=0; fi<xferFields; ++fi) {
				Domain_member src = fieldData[fi] ;
				for (Index_t i=0; i<dx; ++i) {
					destAddr[i] = (domain.*src)(dx*(dy-1) + dx*dy*(dz-1) + i) ;
				}
				destAddr += dx ;
			}
			destAddr -= xferFields*dx ;
			addr[ai++] = try_MPI_Isend(destAddr, xferFields*dx, baseType, toRank, msgType,
					MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg]) ;
			++emsg ;
		}

		if (colMax && planeMax && doSend) {
			int toRank = myRank + domain.tp()*domain.tp() + 1 ;
			destAddr = &domain.commDataSend[pmsg * maxPlaneComm +
				emsg * maxEdgeComm] ;
			for (Index_t fi=0; fi<xferFields; ++fi) {
				Domain_member src = fieldData[fi] ;
				for (Index_t i=0; i<dy; ++i) {
					destAddr[i] = (domain.*src)(dx*dy*(dz-1) + dx - 1 + i*dx) ;
				}
				destAddr += dy ;
			}
			destAddr -= xferFields*dy ;
			addr[ai++] = try_MPI_Isend(destAddr, xferFields*dy, baseType, toRank, msgType,
					MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg]) ;
			++emsg ;
		}

		if (rowMax && colMin && doSend) {
			int toRank = myRank + domain.tp() - 1 ;
			destAddr = &domain.commDataSend[pmsg * maxPlaneComm +
				emsg * maxEdgeComm] ;
			for (Index_t fi=0; fi<xferFields; ++fi) {
				Domain_member src = fieldData[fi] ;
				for (Index_t i=0; i<dz; ++i) {
					destAddr[i] = (domain.*src)(dx*(dy-1) + i*dx*dy) ;
				}
				destAddr += dz ;
			}
			destAddr -= xferFields*dz ;
			addr[ai++] = try_MPI_Isend(destAddr, xferFields*dz, baseType, toRank, msgType,
					MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg]) ;
			++emsg ;
		}

		if (rowMin && planeMax && doSend) {
			int toRank = myRank + domain.tp()*domain.tp() - domain.tp() ;
			destAddr = &domain.commDataSend[pmsg * maxPlaneComm +
				emsg * maxEdgeComm] ;
			for (Index_t fi=0; fi<xferFields; ++fi) {
				Domain_member src = fieldData[fi] ;
				for (Index_t i=0; i<dx; ++i) {
					destAddr[i] = (domain.*src)(dx*dy*(dz-1) + i) ;
				}
				destAddr += dx ;
			}
			destAddr -= xferFields*dx ;
			addr[ai++] = try_MPI_Isend(destAddr, xferFields*dx, baseType, toRank, msgType,
					MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg]) ;
			++emsg ;
		}

		if (colMin && planeMax && doSend) {
			int toRank = myRank + domain.tp()*domain.tp() - 1 ;
			destAddr = &domain.commDataSend[pmsg * maxPlaneComm +
				emsg * maxEdgeComm] ;
			for (Index_t fi=0; fi<xferFields; ++fi) {
				Domain_member src = fieldData[fi] ;
				for (Index_t i=0; i<dy; ++i) {
					destAddr[i] = (domain.*src)(dx*dy*(dz-1) + i*dx) ;
				}
				destAddr += dy ;
			}
			destAddr -= xferFields*dy ;
			addr[ai++] = try_MPI_Isend(destAddr, xferFields*dy, baseType, toRank, msgType,
					MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg]) ;
			++emsg ;
		}

		if (rowMin && colMax) {
			int toRank = myRank - domain.tp() + 1 ;
			destAddr = &domain.commDataSend[pmsg * maxPlaneComm +
				emsg * maxEdgeComm] ;
			for (Index_t fi=0; fi<xferFields; ++fi) {
				Domain_member src = fieldData[fi] ;
				for (Index_t i=0; i<dz; ++i) {
					destAddr[i] = (domain.*src)(dx - 1 + i*dx*dy) ;
				}
				destAddr += dz ;
			}
			destAddr -= xferFields*dz ;
			addr[ai++] = try_MPI_Isend(destAddr, xferFields*dz, baseType, toRank, msgType,
					MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg]) ;
			++emsg ;
		}

		if (rowMax && planeMin) {
			int toRank = myRank - domain.tp()*domain.tp() + domain.tp() ;
			destAddr = &domain.commDataSend[pmsg * maxPlaneComm +
				emsg * maxEdgeComm] ;
			for (Index_t fi=0; fi<xferFields; ++fi) {
				Domain_member src = fieldData[fi] ;
				for (Index_t i=0; i<dx; ++i) {
					destAddr[i] = (domain.*src)(dx*(dy - 1) + i) ;
				}
				destAddr += dx ;
			}
			destAddr -= xferFields*dx ;
			addr[ai++] = try_MPI_Isend(destAddr, xferFields*dx, baseType, toRank, msgType,
					MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg]) ;
			++emsg ;
		}

		if (colMax && planeMin) {
			int toRank = myRank - domain.tp()*domain.tp() + 1 ;
			destAddr = &domain.commDataSend[pmsg * maxPlaneComm +
				emsg * maxEdgeComm] ;
			for (Index_t fi=0; fi<xferFields; ++fi) {
				Domain_member src = fieldData[fi] ;
				for (Index_t i=0; i<dy; ++i) {
					destAddr[i] = (domain.*src)(dx - 1 + i*dx) ;
				}
				destAddr += dy ;
			}
			destAddr -= xferFields*dy ;
			addr[ai++] = try_MPI_Isend(destAddr, xferFields*dy, baseType, toRank, msgType,
					MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg]) ;
			++emsg ;
		}

		if (rowMin && colMin && planeMin) {
			/* corner at domain logical coord (0, 0, 0) */
			int toRank = myRank - domain.tp()*domain.tp() - domain.tp() - 1 ;
			Real_t *comBuf = &domain.commDataSend[pmsg * maxPlaneComm +
				emsg * maxEdgeComm +
				cmsg * CACHE_COHERENCE_PAD_REAL] ;
			for (Index_t fi=0; fi<xferFields; ++fi) {
				comBuf[fi] = (domain.*fieldData[fi])(0) ;
			}
			addr[ai++] = try_MPI_Isend(comBuf, xferFields, baseType, toRank, msgType,
					MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg+cmsg]) ;
			++cmsg ;
		}
		if (rowMin && colMin && planeMax && doSend) {
			/* corner at domain logical coord (0, 0, 1) */
			int toRank = myRank + domain.tp()*domain.tp() - domain.tp() - 1 ;
			Real_t *comBuf = &domain.commDataSend[pmsg * maxPlaneComm +
				emsg * maxEdgeComm +
				cmsg * CACHE_COHERENCE_PAD_REAL] ;
			Index_t idx = dx*dy*(dz - 1) ;
			for (Index_t fi=0; fi<xferFields; ++fi) {
				comBuf[fi] = (domain.*fieldData[fi])(idx) ;
			}
			addr[ai++] = try_MPI_Isend(comBuf, xferFields, baseType, toRank, msgType,
					MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg+cmsg]) ;
			++cmsg ;
		}
		if (rowMin && colMax && planeMin) {
			/* corner at domain logical coord (1, 0, 0) */
			int toRank = myRank - domain.tp()*domain.tp() - domain.tp() + 1 ;
			Real_t *comBuf = &domain.commDataSend[pmsg * maxPlaneComm +
				emsg * maxEdgeComm +
				cmsg * CACHE_COHERENCE_PAD_REAL] ;
			Index_t idx = dx - 1 ;
			for (Index_t fi=0; fi<xferFields; ++fi) {
				comBuf[fi] = (domain.*fieldData[fi])(idx) ;
			}
			addr[ai++] = try_MPI_Isend(comBuf, xferFields, baseType, toRank, msgType,
					MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg+cmsg]) ;
			++cmsg ;
		}
		if (rowMin && colMax && planeMax && doSend) {
			/* corner at domain logical coord (1, 0, 1) */
			int toRank = myRank + domain.tp()*domain.tp() - domain.tp() + 1 ;
			Real_t *comBuf = &domain.commDataSend[pmsg * maxPlaneComm +
				emsg * maxEdgeComm +
				cmsg * CACHE_COHERENCE_PAD_REAL] ;
			Index_t idx = dx*dy*(dz - 1) + (dx - 1) ;
			for (Index_t fi=0; fi<xferFields; ++fi) {
				comBuf[fi] = (domain.*fieldData[fi])(idx) ;
			}
			addr[ai++] = try_MPI_Isend(comBuf, xferFields, baseType, toRank, msgType,
					MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg+cmsg]) ;
			++cmsg ;
		}
		if (rowMax && colMin && planeMin) {
			/* corner at domain logical coord (0, 1, 0) */
			int toRank = myRank - domain.tp()*domain.tp() + domain.tp() - 1 ;
			Real_t *comBuf = &domain.commDataSend[pmsg * maxPlaneComm +
				emsg * maxEdgeComm +
				cmsg * CACHE_COHERENCE_PAD_REAL] ;
			Index_t idx = dx*(dy - 1) ;
			for (Index_t fi=0; fi<xferFields; ++fi) {
				comBuf[fi] = (domain.*fieldData[fi])(idx) ;
			}
			addr[ai++] = try_MPI_Isend(comBuf, xferFields, baseType, toRank, msgType,
					MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg+cmsg]) ;
			++cmsg ;
		}
		if (rowMax && colMin && planeMax && doSend) {
			/* corner at domain logical coord (0, 1, 1) */
			int toRank = myRank + domain.tp()*domain.tp() + domain.tp() - 1 ;
			Real_t *comBuf = &domain.commDataSend[pmsg * maxPlaneComm +
				emsg * maxEdgeComm +
				cmsg * CACHE_COHERENCE_PAD_REAL] ;
			Index_t idx = dx*dy*(dz - 1) + dx*(dy - 1) ;
			for (Index_t fi=0; fi<xferFields; ++fi) {
				comBuf[fi] = (domain.*fieldData[fi])(idx) ;
			}
			addr[ai++] = try_MPI_Isend(comBuf, xferFields, baseType, toRank, msgType,
					MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg+cmsg]) ;
			++cmsg ;
		}
		if (rowMax && colMax && planeMin) {
			/* corner at domain logical coord (1, 1, 0) */
			int toRank = myRank - domain.tp()*domain.tp() + domain.tp() + 1 ;
			Real_t *comBuf = &domain.commDataSend[pmsg * maxPlaneComm +
				emsg * maxEdgeComm +
				cmsg * CACHE_COHERENCE_PAD_REAL] ;
			Index_t idx = dx*dy - 1 ;
			for (Index_t fi=0; fi<xferFields; ++fi) {
				comBuf[fi] = (domain.*fieldData[fi])(idx) ;
			}
			addr[ai++] = try_MPI_Isend(comBuf, xferFields, baseType, toRank, msgType,
					MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg+cmsg]) ;
			++cmsg ;
		}
		if (rowMax && colMax && planeMax && doSend) {
			/* corner at domain logical coord (1, 1, 1) */
			int toRank = myRank + domain.tp()*domain.tp() + domain.tp() + 1 ;
			Real_t *comBuf = &domain.commDataSend[pmsg * maxPlaneComm +
				emsg * maxEdgeComm +
				cmsg * CACHE_COHERENCE_PAD_REAL] ;
			Index_t idx = dx*dy*dz - 1 ;
			for (Index_t fi=0; fi<xferFields; ++fi) {
				comBuf[fi] = (domain.*fieldData[fi])(idx) ;
			}
			addr[ai++] = try_MPI_Isend(comBuf, xferFields, baseType, toRank, msgType,
					MPI_COMM_WORLD, &domain.sendRequest[pmsg+emsg+cmsg]) ;
			++cmsg ;
		}
	}

	MPI_Waitall(26, domain.sendRequest, status);

	for( int i = 0; i < 100; i++ ){
		if( addr[i] != NULL )
			free(addr[i]);
	}

	free(addr);

}

/******************************************/

void CommSBN(Domain& domain, Int_t xferFields, Domain_member *fieldData) {

	if (domain.numRanks() == 1)
		return ;

	/* summation order should be from smallest value to largest */
	/* or we could try out kahan summation! */

	int myRank ;
	Index_t maxPlaneComm = xferFields * domain.maxPlaneSize() ;
	Index_t maxEdgeComm  = xferFields * domain.maxEdgeSize() ;
	Index_t pmsg = 0 ; /* plane comm msg */
	Index_t emsg = 0 ; /* edge comm msg */
	Index_t cmsg = 0 ; /* corner comm msg */
	Index_t dx = domain.sizeX() + 1 ;
	Index_t dy = domain.sizeY() + 1 ;
	Index_t dz = domain.sizeZ() + 1 ;
	MPI_Status status ;
	Real_t *srcAddr ;
	Index_t rowMin, rowMax, colMin, colMax, planeMin, planeMax ;
	/* assume communication to 6 neighbors by default */
	rowMin = rowMax = colMin = colMax = planeMin = planeMax = 1 ;
	if (domain.rowLoc() == 0) {
		rowMin = 0 ;
	}
	if (domain.rowLoc() == (domain.tp()-1)) {
		rowMax = 0 ;
	}
	if (domain.colLoc() == 0) {
		colMin = 0 ;
	}
	if (domain.colLoc() == (domain.tp()-1)) {
		colMax = 0 ;
	}
	if (domain.planeLoc() == 0) {
		planeMin = 0 ;
	}
	if (domain.planeLoc() == (domain.tp()-1)) {
		planeMax = 0 ;
	}

	MPI_Comm_rank(MPI_COMM_WORLD, &myRank) ;

	if (planeMin | planeMax) {
		/* ASSUMING ONE DOMAIN PER RANK, CONSTANT BLOCK SIZE HERE */
		Index_t opCount = dx * dy ;

		if (planeMin) {
			/* contiguous memory */
			srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm] ;
			MPI_Wait(&domain.recvRequest[pmsg], &status) ;

			for (Index_t fi=0 ; fi<xferFields; ++fi) {
				Domain_member dest = fieldData[fi] ;
				for (Index_t i=0; i<opCount; ++i) {
					(domain.*dest)(i) += srcAddr[i] ;
				}
				srcAddr += opCount ;
			}
			++pmsg ;
		}
		if (planeMax) {
			/* contiguous memory */
			srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm] ;
			MPI_Wait(&domain.recvRequest[pmsg], &status) ;
			for (Index_t fi=0 ; fi<xferFields; ++fi) {
				Domain_member dest = fieldData[fi] ;
				for (Index_t i=0; i<opCount; ++i) {
					(domain.*dest)(dx*dy*(dz - 1) + i) += srcAddr[i] ;
				}
				srcAddr += opCount ;
			}
			++pmsg ;
		}
	}

	if (rowMin | rowMax) {
		/* ASSUMING ONE DOMAIN PER RANK, CONSTANT BLOCK SIZE HERE */
		Index_t opCount = dx * dz ;

		if (rowMin) {
			/* contiguous memory */
			srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm] ;
			MPI_Wait(&domain.recvRequest[pmsg], &status) ;
			for (Index_t fi=0 ; fi<xferFields; ++fi) {
				Domain_member dest = fieldData[fi] ;
				for (Index_t i=0; i<dz; ++i) {
					for (Index_t j=0; j<dx; ++j) {
						(domain.*dest)(i*dx*dy + j) += srcAddr[i*dx + j] ;
					}
				}
				srcAddr += opCount ;
			}
			++pmsg ;
		}
		if (rowMax) {
			/* contiguous memory */
			srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm] ;
			MPI_Wait(&domain.recvRequest[pmsg], &status) ;
			for (Index_t fi=0 ; fi<xferFields; ++fi) {
				Domain_member dest = fieldData[fi] ;
				for (Index_t i=0; i<dz; ++i) {
					for (Index_t j=0; j<dx; ++j) {
						(domain.*dest)(dx*(dy - 1) + i*dx*dy + j) += srcAddr[i*dx + j] ;
					}
				}
				srcAddr += opCount ;
			}
			++pmsg ;
		}
	}
	if (colMin | colMax) {
		/* ASSUMING ONE DOMAIN PER RANK, CONSTANT BLOCK SIZE HERE */
		Index_t opCount = dy * dz ;

		if (colMin) {
			/* contiguous memory */
			srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm] ;
			MPI_Wait(&domain.recvRequest[pmsg], &status) ;
			for (Index_t fi=0 ; fi<xferFields; ++fi) {
				Domain_member dest = fieldData[fi] ;
				for (Index_t i=0; i<dz; ++i) {
					for (Index_t j=0; j<dy; ++j) {
						(domain.*dest)(i*dx*dy + j*dx) += srcAddr[i*dy + j] ;
					}
				}
				srcAddr += opCount ;
			}
			++pmsg ;
		}
		if (colMax) {
			/* contiguous memory */
			srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm] ;
			MPI_Wait(&domain.recvRequest[pmsg], &status) ;
			for (Index_t fi=0 ; fi<xferFields; ++fi) {
				Domain_member dest = fieldData[fi] ;
				for (Index_t i=0; i<dz; ++i) {
					for (Index_t j=0; j<dy; ++j) {
						(domain.*dest)(dx - 1 + i*dx*dy + j*dx) += srcAddr[i*dy + j] ;
					}
				}
				srcAddr += opCount ;
			}
			++pmsg ;
		}
	}

	if (rowMin & colMin) {
		srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm +
			emsg * maxEdgeComm] ;
		MPI_Wait(&domain.recvRequest[pmsg+emsg], &status) ;
		for (Index_t fi=0 ; fi<xferFields; ++fi) {
			Domain_member dest = fieldData[fi] ;
			for (Index_t i=0; i<dz; ++i) {
				(domain.*dest)(i*dx*dy) += srcAddr[i] ;
			}
			srcAddr += dz ;
		}
		++emsg ;
	}

	if (rowMin & planeMin) {
		srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm +
			emsg * maxEdgeComm] ;
		MPI_Wait(&domain.recvRequest[pmsg+emsg], &status) ;
		for (Index_t fi=0 ; fi<xferFields; ++fi) {
			Domain_member dest = fieldData[fi] ;
			for (Index_t i=0; i<dx; ++i) {
				(domain.*dest)(i) += srcAddr[i] ;
			}
			srcAddr += dx ;
		}
		++emsg ;
	}

	if (colMin & planeMin) {
		srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm +
			emsg * maxEdgeComm] ;
		MPI_Wait(&domain.recvRequest[pmsg+emsg], &status) ;
		for (Index_t fi=0 ; fi<xferFields; ++fi) {
			Domain_member dest = fieldData[fi] ;
			for (Index_t i=0; i<dy; ++i) {
				(domain.*dest)(i*dx) += srcAddr[i] ;
			}
			srcAddr += dy ;
		}
		++emsg ;
	}

	if (rowMax & colMax) {
		srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm +
			emsg * maxEdgeComm] ;
		MPI_Wait(&domain.recvRequest[pmsg+emsg], &status) ;
		for (Index_t fi=0 ; fi<xferFields; ++fi) {
			Domain_member dest = fieldData[fi] ;
			for (Index_t i=0; i<dz; ++i) {
				(domain.*dest)(dx*dy - 1 + i*dx*dy) += srcAddr[i] ;
			}
			srcAddr += dz ;
		}
		++emsg ;
	}

	if (rowMax & planeMax) {
		srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm +
			emsg * maxEdgeComm] ;
		MPI_Wait(&domain.recvRequest[pmsg+emsg], &status) ;
		for (Index_t fi=0 ; fi<xferFields; ++fi) {
			Domain_member dest = fieldData[fi] ;
			for (Index_t i=0; i<dx; ++i) {
				(domain.*dest)(dx*(dy-1) + dx*dy*(dz-1) + i) += srcAddr[i] ;
			}
			srcAddr += dx ;
		}
		++emsg ;
	}

	if (colMax & planeMax) {
		srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm +
			emsg * maxEdgeComm] ;
		MPI_Wait(&domain.recvRequest[pmsg+emsg], &status) ;
		for (Index_t fi=0 ; fi<xferFields; ++fi) {
			Domain_member dest = fieldData[fi] ;
			for (Index_t i=0; i<dy; ++i) {
				(domain.*dest)(dx*dy*(dz-1) + dx - 1 + i*dx) += srcAddr[i] ;
			}
			srcAddr += dy ;
		}
		++emsg ;
	}

	if (rowMax & colMin) {
		srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm +
			emsg * maxEdgeComm] ;
		MPI_Wait(&domain.recvRequest[pmsg+emsg], &status) ;
		for (Index_t fi=0 ; fi<xferFields; ++fi) {
			Domain_member dest = fieldData[fi] ;
			for (Index_t i=0; i<dz; ++i) {
				(domain.*dest)(dx*(dy-1) + i*dx*dy) += srcAddr[i] ;
			}
			srcAddr += dz ;
		}
		++emsg ;
	}

	if (rowMin & planeMax) {
		srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm +
			emsg * maxEdgeComm] ;
		MPI_Wait(&domain.recvRequest[pmsg+emsg], &status) ;
		for (Index_t fi=0 ; fi<xferFields; ++fi) {
			Domain_member dest = fieldData[fi] ;
			for (Index_t i=0; i<dx; ++i) {
				(domain.*dest)(dx*dy*(dz-1) + i) += srcAddr[i] ;
			}
			srcAddr += dx ;
		}
		++emsg ;
	}

	if (colMin & planeMax) {
		srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm +
			emsg * maxEdgeComm] ;
		MPI_Wait(&domain.recvRequest[pmsg+emsg], &status) ;
		for (Index_t fi=0 ; fi<xferFields; ++fi) {
			Domain_member dest = fieldData[fi] ;
			for (Index_t i=0; i<dy; ++i) {
				(domain.*dest)(dx*dy*(dz-1) + i*dx) += srcAddr[i] ;
			}
			srcAddr += dy ;
		}
		++emsg ;
	}

	if (rowMin & colMax) {
		srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm +
			emsg * maxEdgeComm] ;
		MPI_Wait(&domain.recvRequest[pmsg+emsg], &status) ;
		for (Index_t fi=0 ; fi<xferFields; ++fi) {
			Domain_member dest = fieldData[fi] ;
			for (Index_t i=0; i<dz; ++i) {
				(domain.*dest)(dx - 1 + i*dx*dy) += srcAddr[i] ;
			}
			srcAddr += dz ;
		}
		++emsg ;
	}

	if (rowMax & planeMin) {
		srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm +
			emsg * maxEdgeComm] ;
		MPI_Wait(&domain.recvRequest[pmsg+emsg], &status) ;
		for (Index_t fi=0 ; fi<xferFields; ++fi) {
			Domain_member dest = fieldData[fi] ;
			for (Index_t i=0; i<dx; ++i) {
				(domain.*dest)(dx*(dy - 1) + i) += srcAddr[i] ;
			}
			srcAddr += dx ;
		}
		++emsg ;
	}

	if (colMax & planeMin) {
		srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm +
			emsg * maxEdgeComm] ;
		MPI_Wait(&domain.recvRequest[pmsg+emsg], &status) ;
		for (Index_t fi=0 ; fi<xferFields; ++fi) {
			Domain_member dest = fieldData[fi] ;
			for (Index_t i=0; i<dy; ++i) {
				(domain.*dest)(dx - 1 + i*dx) += srcAddr[i] ;
			}
			srcAddr += dy ;
		}
		++emsg ;
	}

	if (rowMin & colMin & planeMin) {
		/* corner at domain logical coord (0, 0, 0) */
		Real_t *comBuf = &domain.commDataRecv[pmsg * maxPlaneComm +
			emsg * maxEdgeComm +
			cmsg * CACHE_COHERENCE_PAD_REAL] ;
		MPI_Wait(&domain.recvRequest[pmsg+emsg+cmsg], &status) ;
		for (Index_t fi=0; fi<xferFields; ++fi) {
			(domain.*fieldData[fi])(0) += comBuf[fi] ;
		}
		++cmsg ;
	}
	if (rowMin & colMin & planeMax) {
		/* corner at domain logical coord (0, 0, 1) */
		Real_t *comBuf = &domain.commDataRecv[pmsg * maxPlaneComm +
			emsg * maxEdgeComm +
			cmsg * CACHE_COHERENCE_PAD_REAL] ;
		Index_t idx = dx*dy*(dz - 1) ;
		MPI_Wait(&domain.recvRequest[pmsg+emsg+cmsg], &status) ;
		for (Index_t fi=0; fi<xferFields; ++fi) {
			(domain.*fieldData[fi])(idx) += comBuf[fi] ;
		}
		++cmsg ;
	}
	if (rowMin & colMax & planeMin) {
		/* corner at domain logical coord (1, 0, 0) */
		Real_t *comBuf = &domain.commDataRecv[pmsg * maxPlaneComm +
			emsg * maxEdgeComm +
			cmsg * CACHE_COHERENCE_PAD_REAL] ;
		Index_t idx = dx - 1 ;
		MPI_Wait(&domain.recvRequest[pmsg+emsg+cmsg], &status) ;
		for (Index_t fi=0; fi<xferFields; ++fi) {
			(domain.*fieldData[fi])(idx) += comBuf[fi] ;
		}
		++cmsg ;
	}
	if (rowMin & colMax & planeMax) {
		/* corner at domain logical coord (1, 0, 1) */
		Real_t *comBuf = &domain.commDataRecv[pmsg * maxPlaneComm +
			emsg * maxEdgeComm +
			cmsg * CACHE_COHERENCE_PAD_REAL] ;
		Index_t idx = dx*dy*(dz - 1) + (dx - 1) ;
		MPI_Wait(&domain.recvRequest[pmsg+emsg+cmsg], &status) ;
		for (Index_t fi=0; fi<xferFields; ++fi) {
			(domain.*fieldData[fi])(idx) += comBuf[fi] ;
		}
		++cmsg ;
	}
	if (rowMax & colMin & planeMin) {
		/* corner at domain logical coord (0, 1, 0) */
		Real_t *comBuf = &domain.commDataRecv[pmsg * maxPlaneComm +
			emsg * maxEdgeComm +
			cmsg * CACHE_COHERENCE_PAD_REAL] ;
		Index_t idx = dx*(dy - 1) ;
		MPI_Wait(&domain.recvRequest[pmsg+emsg+cmsg], &status) ;
		for (Index_t fi=0; fi<xferFields; ++fi) {
			(domain.*fieldData[fi])(idx) += comBuf[fi] ;
		}
		++cmsg ;
	}
	if (rowMax & colMin & planeMax) {
		/* corner at domain logical coord (0, 1, 1) */
		Real_t *comBuf = &domain.commDataRecv[pmsg * maxPlaneComm +
			emsg * maxEdgeComm +
			cmsg * CACHE_COHERENCE_PAD_REAL] ;
		Index_t idx = dx*dy*(dz - 1) + dx*(dy - 1) ;
		MPI_Wait(&domain.recvRequest[pmsg+emsg+cmsg], &status) ;
		for (Index_t fi=0; fi<xferFields; ++fi) {
			(domain.*fieldData[fi])(idx) += comBuf[fi] ;
		}
		++cmsg ;
	}
	if (rowMax & colMax & planeMin) {
		/* corner at domain logical coord (1, 1, 0) */
		Real_t *comBuf = &domain.commDataRecv[pmsg * maxPlaneComm +
			emsg * maxEdgeComm +
			cmsg * CACHE_COHERENCE_PAD_REAL] ;
		Index_t idx = dx*dy - 1 ;
		MPI_Wait(&domain.recvRequest[pmsg+emsg+cmsg], &status) ;
		for (Index_t fi=0; fi<xferFields; ++fi) {
			(domain.*fieldData[fi])(idx) += comBuf[fi] ;
		}
		++cmsg ;
	}
	if (rowMax & colMax & planeMax) {
		/* corner at domain logical coord (1, 1, 1) */
		Real_t *comBuf = &domain.commDataRecv[pmsg * maxPlaneComm +
			emsg * maxEdgeComm +
			cmsg * CACHE_COHERENCE_PAD_REAL] ;
		Index_t idx = dx*dy*dz - 1 ;
		MPI_Wait(&domain.recvRequest[pmsg+emsg+cmsg], &status) ;
		for (Index_t fi=0; fi<xferFields; ++fi) {
			(domain.*fieldData[fi])(idx) += comBuf[fi] ;
		}
		++cmsg ;
	}
}

/******************************************/

void CommSyncPosVel(Domain& domain) {

	if (domain.numRanks() == 1)
		return ;

	int myRank ;
	bool doRecv = false ;
	Index_t xferFields = 6 ; /* x, y, z, xd, yd, zd */
	Domain_member fieldData[6] ;
	Index_t maxPlaneComm = xferFields * domain.maxPlaneSize() ;
	Index_t maxEdgeComm  = xferFields * domain.maxEdgeSize() ;
	Index_t pmsg = 0 ; /* plane comm msg */
	Index_t emsg = 0 ; /* edge comm msg */
	Index_t cmsg = 0 ; /* corner comm msg */
	Index_t dx = domain.sizeX() + 1 ;
	Index_t dy = domain.sizeY() + 1 ;
	Index_t dz = domain.sizeZ() + 1 ;
	MPI_Status status ;
	Real_t *srcAddr ;
	bool rowMin, rowMax, colMin, colMax, planeMin, planeMax ;

	/* assume communication to 6 neighbors by default */
	rowMin = rowMax = colMin = colMax = planeMin = planeMax = true ;
	if (domain.rowLoc() == 0) {
		rowMin = false ;
	}
	if (domain.rowLoc() == (domain.tp()-1)) {
		rowMax = false ;
	}
	if (domain.colLoc() == 0) {
		colMin = false ;
	}
	if (domain.colLoc() == (domain.tp()-1)) {
		colMax = false ;
	}
	if (domain.planeLoc() == 0) {
		planeMin = false ;
	}
	if (domain.planeLoc() == (domain.tp()-1)) {
		planeMax = false ;
	}

	fieldData[0] = &Domain::x ;
	fieldData[1] = &Domain::y ;
	fieldData[2] = &Domain::z ;
	fieldData[3] = &Domain::xd ;
	fieldData[4] = &Domain::yd ;
	fieldData[5] = &Domain::zd ;

	MPI_Comm_rank(MPI_COMM_WORLD, &myRank) ;

	if (planeMin | planeMax) {
		/* ASSUMING ONE DOMAIN PER RANK, CONSTANT BLOCK SIZE HERE */
		Index_t opCount = dx * dy ;

		if (planeMin && doRecv) {
			/* contiguous memory */
			srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm] ;
			MPI_Wait(&domain.recvRequest[pmsg], &status) ;
			for (Index_t fi=0 ; fi<xferFields; ++fi) {
				Domain_member dest = fieldData[fi] ;
				for (Index_t i=0; i<opCount; ++i) {
					(domain.*dest)(i) = srcAddr[i] ;
				}
				srcAddr += opCount ;
			}
			++pmsg ;
		}
		if (planeMax) {
			/* contiguous memory */
			srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm] ;
			MPI_Wait(&domain.recvRequest[pmsg], &status) ;
			for (Index_t fi=0 ; fi<xferFields; ++fi) {
				Domain_member dest = fieldData[fi] ;
				for (Index_t i=0; i<opCount; ++i) {
					(domain.*dest)(dx*dy*(dz - 1) + i) = srcAddr[i] ;
				}
				srcAddr += opCount ;
			}
			++pmsg ;
		}
	}

	if (rowMin | rowMax) {
		/* ASSUMING ONE DOMAIN PER RANK, CONSTANT BLOCK SIZE HERE */
		Index_t opCount = dx * dz ;

		if (rowMin && doRecv) {
			/* contiguous memory */
			srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm] ;
			MPI_Wait(&domain.recvRequest[pmsg], &status) ;
			for (Index_t fi=0 ; fi<xferFields; ++fi) {
				Domain_member dest = fieldData[fi] ;
				for (Index_t i=0; i<dz; ++i) {
					for (Index_t j=0; j<dx; ++j) {
						(domain.*dest)(i*dx*dy + j) = srcAddr[i*dx + j] ;
					}
				}
				srcAddr += opCount ;
			}
			++pmsg ;
		}
		if (rowMax) {
			/* contiguous memory */
			srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm] ;
			MPI_Wait(&domain.recvRequest[pmsg], &status) ;
			for (Index_t fi=0 ; fi<xferFields; ++fi) {
				Domain_member dest = fieldData[fi] ;
				for (Index_t i=0; i<dz; ++i) {
					for (Index_t j=0; j<dx; ++j) {
						(domain.*dest)(dx*(dy - 1) + i*dx*dy + j) = srcAddr[i*dx + j] ;
					}
				}
				srcAddr += opCount ;
			}
			++pmsg ;
		}
	}

	if (colMin | colMax) {
		/* ASSUMING ONE DOMAIN PER RANK, CONSTANT BLOCK SIZE HERE */
		Index_t opCount = dy * dz ;

		if (colMin && doRecv) {
			/* contiguous memory */
			srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm] ;
			MPI_Wait(&domain.recvRequest[pmsg], &status) ;
			for (Index_t fi=0 ; fi<xferFields; ++fi) {
				Domain_member dest = fieldData[fi] ;
				for (Index_t i=0; i<dz; ++i) {
					for (Index_t j=0; j<dy; ++j) {
						(domain.*dest)(i*dx*dy + j*dx) = srcAddr[i*dy + j] ;
					}
				}
				srcAddr += opCount ;
			}
			++pmsg ;
		}
		if (colMax) {
			/* contiguous memory */
			srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm] ;
			MPI_Wait(&domain.recvRequest[pmsg], &status) ;
			for (Index_t fi=0 ; fi<xferFields; ++fi) {
				Domain_member dest = fieldData[fi] ;
				for (Index_t i=0; i<dz; ++i) {
					for (Index_t j=0; j<dy; ++j) {
						(domain.*dest)(dx - 1 + i*dx*dy + j*dx) = srcAddr[i*dy + j] ;
					}
				}
				srcAddr += opCount ;
			}
			++pmsg ;
		}
	}

	if (rowMin && colMin && doRecv) {
		srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm +
			emsg * maxEdgeComm] ;
		MPI_Wait(&domain.recvRequest[pmsg+emsg], &status) ;
		for (Index_t fi=0 ; fi<xferFields; ++fi) {
			Domain_member dest = fieldData[fi] ;
			for (Index_t i=0; i<dz; ++i) {
				(domain.*dest)(i*dx*dy) = srcAddr[i] ;
			}
			srcAddr += dz ;
		}
		++emsg ;
	}

	if (rowMin && planeMin && doRecv) {
		srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm +
			emsg * maxEdgeComm] ;
		MPI_Wait(&domain.recvRequest[pmsg+emsg], &status) ;
		for (Index_t fi=0 ; fi<xferFields; ++fi) {
			Domain_member dest = fieldData[fi] ;
			for (Index_t i=0; i<dx; ++i) {
				(domain.*dest)(i) = srcAddr[i] ;
			}
			srcAddr += dx ;
		}
		++emsg ;
	}

	if (colMin && planeMin && doRecv) {
		srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm +
			emsg * maxEdgeComm] ;
		MPI_Wait(&domain.recvRequest[pmsg+emsg], &status) ;
		for (Index_t fi=0 ; fi<xferFields; ++fi) {
			Domain_member dest = fieldData[fi] ;
			for (Index_t i=0; i<dy; ++i) {
				(domain.*dest)(i*dx) = srcAddr[i] ;
			}
			srcAddr += dy ;
		}
		++emsg ;
	}

	if (rowMax && colMax) {
		srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm +
			emsg * maxEdgeComm] ;
		MPI_Wait(&domain.recvRequest[pmsg+emsg], &status) ;
		for (Index_t fi=0 ; fi<xferFields; ++fi) {
			Domain_member dest = fieldData[fi] ;
			for (Index_t i=0; i<dz; ++i) {
				(domain.*dest)(dx*dy - 1 + i*dx*dy) = srcAddr[i] ;
			}
			srcAddr += dz ;
		}
		++emsg ;
	}

	if (rowMax && planeMax) {
		srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm +
			emsg * maxEdgeComm] ;
		MPI_Wait(&domain.recvRequest[pmsg+emsg], &status) ;
		for (Index_t fi=0 ; fi<xferFields; ++fi) {
			Domain_member dest = fieldData[fi] ;
			for (Index_t i=0; i<dx; ++i) {
				(domain.*dest)(dx*(dy-1) + dx*dy*(dz-1) + i) = srcAddr[i] ;
			}
			srcAddr += dx ;
		}
		++emsg ;
	}

	if (colMax && planeMax) {
		srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm +
			emsg * maxEdgeComm] ;
		MPI_Wait(&domain.recvRequest[pmsg+emsg], &status) ;
		for (Index_t fi=0 ; fi<xferFields; ++fi) {
			Domain_member dest = fieldData[fi] ;
			for (Index_t i=0; i<dy; ++i) {
				(domain.*dest)(dx*dy*(dz-1) + dx - 1 + i*dx) = srcAddr[i] ;
			}
			srcAddr += dy ;
		}
		++emsg ;
	}

	if (rowMax && colMin) {
		srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm +
			emsg * maxEdgeComm] ;
		MPI_Wait(&domain.recvRequest[pmsg+emsg], &status) ;
		for (Index_t fi=0 ; fi<xferFields; ++fi) {
			Domain_member dest = fieldData[fi] ;
			for (Index_t i=0; i<dz; ++i) {
				(domain.*dest)(dx*(dy-1) + i*dx*dy) = srcAddr[i] ;
			}
			srcAddr += dz ;
		}
		++emsg ;
	}

	if (rowMin && planeMax) {
		srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm +
			emsg * maxEdgeComm] ;
		MPI_Wait(&domain.recvRequest[pmsg+emsg], &status) ;
		for (Index_t fi=0 ; fi<xferFields; ++fi) {
			Domain_member dest = fieldData[fi] ;
			for (Index_t i=0; i<dx; ++i) {
				(domain.*dest)(dx*dy*(dz-1) + i) = srcAddr[i] ;
			}
			srcAddr += dx ;
		}
		++emsg ;
	}

	if (colMin && planeMax) {
		srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm +
			emsg * maxEdgeComm] ;
		MPI_Wait(&domain.recvRequest[pmsg+emsg], &status) ;
		for (Index_t fi=0 ; fi<xferFields; ++fi) {
			Domain_member dest = fieldData[fi] ;
			for (Index_t i=0; i<dy; ++i) {
				(domain.*dest)(dx*dy*(dz-1) + i*dx) = srcAddr[i] ;
			}
			srcAddr += dy ;
		}
		++emsg ;
	}

	if (rowMin && colMax && doRecv) {
		srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm +
			emsg * maxEdgeComm] ;
		MPI_Wait(&domain.recvRequest[pmsg+emsg], &status) ;
		for (Index_t fi=0 ; fi<xferFields; ++fi) {
			Domain_member dest = fieldData[fi] ;
			for (Index_t i=0; i<dz; ++i) {
				(domain.*dest)(dx - 1 + i*dx*dy) = srcAddr[i] ;
			}
			srcAddr += dz ;
		}
		++emsg ;
	}

	if (rowMax && planeMin && doRecv) {
		srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm +
			emsg * maxEdgeComm] ;
		MPI_Wait(&domain.recvRequest[pmsg+emsg], &status) ;
		for (Index_t fi=0 ; fi<xferFields; ++fi) {
			Domain_member dest = fieldData[fi] ;
			for (Index_t i=0; i<dx; ++i) {
				(domain.*dest)(dx*(dy - 1) + i) = srcAddr[i] ;
			}
			srcAddr += dx ;
		}
		++emsg ;
	}

	if (colMax && planeMin && doRecv) {
		srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm +
			emsg * maxEdgeComm] ;
		MPI_Wait(&domain.recvRequest[pmsg+emsg], &status) ;
		for (Index_t fi=0 ; fi<xferFields; ++fi) {
			Domain_member dest = fieldData[fi] ;
			for (Index_t i=0; i<dy; ++i) {
				(domain.*dest)(dx - 1 + i*dx) = srcAddr[i] ;
			}
			srcAddr += dy ;
		}
		++emsg ;
	}


	if (rowMin && colMin && planeMin && doRecv) {
		/* corner at domain logical coord (0, 0, 0) */
		Real_t *comBuf = &domain.commDataRecv[pmsg * maxPlaneComm +
			emsg * maxEdgeComm +
			cmsg * CACHE_COHERENCE_PAD_REAL] ;
		MPI_Wait(&domain.recvRequest[pmsg+emsg+cmsg], &status) ;
		for (Index_t fi=0; fi<xferFields; ++fi) {
			(domain.*fieldData[fi])(0) = comBuf[fi] ;
		}
		++cmsg ;
	}
	if (rowMin && colMin && planeMax) {
		/* corner at domain logical coord (0, 0, 1) */
		Real_t *comBuf = &domain.commDataRecv[pmsg * maxPlaneComm +
			emsg * maxEdgeComm +
			cmsg * CACHE_COHERENCE_PAD_REAL] ;
		Index_t idx = dx*dy*(dz - 1) ;
		MPI_Wait(&domain.recvRequest[pmsg+emsg+cmsg], &status) ;
		for (Index_t fi=0; fi<xferFields; ++fi) {
			(domain.*fieldData[fi])(idx) = comBuf[fi] ;
		}
		++cmsg ;
	}
	if (rowMin && colMax && planeMin && doRecv) {
		/* corner at domain logical coord (1, 0, 0) */
		Real_t *comBuf = &domain.commDataRecv[pmsg * maxPlaneComm +
			emsg * maxEdgeComm +
			cmsg * CACHE_COHERENCE_PAD_REAL] ;
		Index_t idx = dx - 1 ;
		MPI_Wait(&domain.recvRequest[pmsg+emsg+cmsg], &status) ;
		for (Index_t fi=0; fi<xferFields; ++fi) {
			(domain.*fieldData[fi])(idx) = comBuf[fi] ;
		}
		++cmsg ;
	}
	if (rowMin && colMax && planeMax) {
		/* corner at domain logical coord (1, 0, 1) */
		Real_t *comBuf = &domain.commDataRecv[pmsg * maxPlaneComm +
			emsg * maxEdgeComm +
			cmsg * CACHE_COHERENCE_PAD_REAL] ;
		Index_t idx = dx*dy*(dz - 1) + (dx - 1) ;
		MPI_Wait(&domain.recvRequest[pmsg+emsg+cmsg], &status) ;
		for (Index_t fi=0; fi<xferFields; ++fi) {
			(domain.*fieldData[fi])(idx) = comBuf[fi] ;
		}
		++cmsg ;
	}
	if (rowMax && colMin && planeMin && doRecv) {
		/* corner at domain logical coord (0, 1, 0) */
		Real_t *comBuf = &domain.commDataRecv[pmsg * maxPlaneComm +
			emsg * maxEdgeComm +
			cmsg * CACHE_COHERENCE_PAD_REAL] ;
		Index_t idx = dx*(dy - 1) ;
		MPI_Wait(&domain.recvRequest[pmsg+emsg+cmsg], &status) ;
		for (Index_t fi=0; fi<xferFields; ++fi) {
			(domain.*fieldData[fi])(idx) = comBuf[fi] ;
		}
		++cmsg ;
	}
	if (rowMax && colMin && planeMax) {
		/* corner at domain logical coord (0, 1, 1) */
		Real_t *comBuf = &domain.commDataRecv[pmsg * maxPlaneComm +
			emsg * maxEdgeComm +
			cmsg * CACHE_COHERENCE_PAD_REAL] ;
		Index_t idx = dx*dy*(dz - 1) + dx*(dy - 1) ;
		MPI_Wait(&domain.recvRequest[pmsg+emsg+cmsg], &status) ;
		for (Index_t fi=0; fi<xferFields; ++fi) {
			(domain.*fieldData[fi])(idx) = comBuf[fi] ;
		}
		++cmsg ;
	}
	if (rowMax && colMax && planeMin && doRecv) {
		/* corner at domain logical coord (1, 1, 0) */
		Real_t *comBuf = &domain.commDataRecv[pmsg * maxPlaneComm +
			emsg * maxEdgeComm +
			cmsg * CACHE_COHERENCE_PAD_REAL] ;
		Index_t idx = dx*dy - 1 ;
		MPI_Wait(&domain.recvRequest[pmsg+emsg+cmsg], &status) ;
		for (Index_t fi=0; fi<xferFields; ++fi) {
			(domain.*fieldData[fi])(idx) = comBuf[fi] ;
		}
		++cmsg ;
	}
	if (rowMax && colMax && planeMax) {
		/* corner at domain logical coord (1, 1, 1) */
		Real_t *comBuf = &domain.commDataRecv[pmsg * maxPlaneComm +
			emsg * maxEdgeComm +
			cmsg * CACHE_COHERENCE_PAD_REAL] ;
		Index_t idx = dx*dy*dz - 1 ;
		MPI_Wait(&domain.recvRequest[pmsg+emsg+cmsg], &status) ;
		for (Index_t fi=0; fi<xferFields; ++fi) {
			(domain.*fieldData[fi])(idx) = comBuf[fi] ;
		}
		++cmsg ;
	}
}

/******************************************/

void CommMonoQ(Domain& domain)
{
	if (domain.numRanks() == 1)
		return ;

	int myRank ;
	Index_t xferFields = 3 ; /* delv_xi, delv_eta, delv_zeta */
	Domain_member fieldData[3] ;
	Index_t fieldOffset[3] ;
	Index_t maxPlaneComm = xferFields * domain.maxPlaneSize() ;
	Index_t pmsg = 0 ; /* plane comm msg */
	Index_t dx = domain.sizeX() ;
	Index_t dy = domain.sizeY() ;
	Index_t dz = domain.sizeZ() ;
	MPI_Status status ;
	Real_t *srcAddr ;
	bool rowMin, rowMax, colMin, colMax, planeMin, planeMax ;
	/* assume communication to 6 neighbors by default */
	rowMin = rowMax = colMin = colMax = planeMin = planeMax = true ;
	if (domain.rowLoc() == 0) {
		rowMin = false ;
	}
	if (domain.rowLoc() == (domain.tp()-1)) {
		rowMax = false ;
	}
	if (domain.colLoc() == 0) {
		colMin = false ;
	}
	if (domain.colLoc() == (domain.tp()-1)) {
		colMax = false ;
	}
	if (domain.planeLoc() == 0) {
		planeMin = false ;
	}
	if (domain.planeLoc() == (domain.tp()-1)) {
		planeMax = false ;
	}

	/* point into ghost data area */
	// fieldData[0] = &(domain.delv_xi(domain.numElem())) ;
	// fieldData[1] = &(domain.delv_eta(domain.numElem())) ;
	// fieldData[2] = &(domain.delv_zeta(domain.numElem())) ;
	fieldData[0] = &Domain::delv_xi ;
	fieldData[1] = &Domain::delv_eta ;
	fieldData[2] = &Domain::delv_zeta ;
	fieldOffset[0] = domain.numElem() ;
	fieldOffset[1] = domain.numElem() ;
	fieldOffset[2] = domain.numElem() ;


	MPI_Comm_rank(MPI_COMM_WORLD, &myRank) ;

	if (planeMin | planeMax) {
		/* ASSUMING ONE DOMAIN PER RANK, CONSTANT BLOCK SIZE HERE */
		Index_t opCount = dx * dy ;

		if (planeMin) {
			/* contiguous memory */
			srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm] ;
			MPI_Wait(&domain.recvRequest[pmsg], &status) ;
			for (Index_t fi=0 ; fi<xferFields; ++fi) {
				Domain_member dest = fieldData[fi] ;
				for (Index_t i=0; i<opCount; ++i) {
					(domain.*dest)(fieldOffset[fi] + i) = srcAddr[i] ;
				}
				srcAddr += opCount ;
				fieldOffset[fi] += opCount ;
			}
			++pmsg ;
		}
		if (planeMax) {
			/* contiguous memory */
			srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm] ;
			MPI_Wait(&domain.recvRequest[pmsg], &status) ;
			for (Index_t fi=0 ; fi<xferFields; ++fi) {
				Domain_member dest = fieldData[fi] ;
				for (Index_t i=0; i<opCount; ++i) {
					(domain.*dest)(fieldOffset[fi] + i) = srcAddr[i] ;
				}
				srcAddr += opCount ;
				fieldOffset[fi] += opCount ;
			}
			++pmsg ;
		}
	}

	if (rowMin | rowMax) {
		/* ASSUMING ONE DOMAIN PER RANK, CONSTANT BLOCK SIZE HERE */
		Index_t opCount = dx * dz ;

		if (rowMin) {
			/* contiguous memory */
			srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm] ;
			MPI_Wait(&domain.recvRequest[pmsg], &status) ;
			for (Index_t fi=0 ; fi<xferFields; ++fi) {
				Domain_member dest = fieldData[fi] ;
				for (Index_t i=0; i<opCount; ++i) {
					(domain.*dest)(fieldOffset[fi] + i) = srcAddr[i] ;
				}
				srcAddr += opCount ;
				fieldOffset[fi] += opCount ;
			}
			++pmsg ;
		}
		if (rowMax) {
			/* contiguous memory */
			srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm] ;
			MPI_Wait(&domain.recvRequest[pmsg], &status) ;
			for (Index_t fi=0 ; fi<xferFields; ++fi) {
				Domain_member dest = fieldData[fi] ;
				for (Index_t i=0; i<opCount; ++i) {
					(domain.*dest)(fieldOffset[fi] + i) = srcAddr[i] ;
				}
				srcAddr += opCount ;
				fieldOffset[fi] += opCount ;
			}
			++pmsg ;
		}
	}
	if (colMin | colMax) {
		/* ASSUMING ONE DOMAIN PER RANK, CONSTANT BLOCK SIZE HERE */
		Index_t opCount = dy * dz ;

		if (colMin) {
			/* contiguous memory */
			srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm] ;
			MPI_Wait(&domain.recvRequest[pmsg], &status) ;
			for (Index_t fi=0 ; fi<xferFields; ++fi) {
				Domain_member dest = fieldData[fi] ;
				for (Index_t i=0; i<opCount; ++i) {
					(domain.*dest)(fieldOffset[fi] + i) = srcAddr[i] ;
				}
				srcAddr += opCount ;
				fieldOffset[fi] += opCount ;
			}
			++pmsg ;
		}
		if (colMax) {
			/* contiguous memory */
			srcAddr = &domain.commDataRecv[pmsg * maxPlaneComm] ;
			MPI_Wait(&domain.recvRequest[pmsg], &status) ;
			for (Index_t fi=0 ; fi<xferFields; ++fi) {
				Domain_member dest = fieldData[fi] ;
				for (Index_t i=0; i<opCount; ++i) {
					(domain.*dest)(fieldOffset[fi] + i) = srcAddr[i] ;
				}
				srcAddr += opCount ;
			}
			++pmsg ;
		}
	}
}

#endif
