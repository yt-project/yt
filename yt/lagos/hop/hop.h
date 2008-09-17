#include "slice.h"
//#define free(A) if(A==NULL)fprintf(stderr,"FREEING DOUBLE\n");fprintf(stderr,"Freeing "#A" ("__FILE__":%d)\n",__LINE__);free(A);

/* ----------------------------------------------------------------------- */
/* The following structures track all the information about the groups */
 
typedef struct groupstruct {
    int npart;          /* Number of particles in the group */
    int npartcum;       /* Cumulative number of particles */
    int nread;          /* Number read so far, also a utility field */
    double compos[3], comvel[3];/* Lists of group CoM position and velocities */
    double comtemp[3];  /* Temporary CoM position */
    int idmerge;        /* The new group ID # after merging.  Not necessarily
                                unique! */
    int rootgroup;	/* The fully traced group id */
} Group;  /* Type Group is defined */
 
typedef struct groupliststruct {
    int npart;          /* Number of particles in the simulation */
    int ngroups;        /* Number of groups in list */
    int nnewgroups;     /* Number of groups after relabeling */
    int npartingroups;  /* Number of particles in groups */
    Group *list;        /* List of groups, zero-offset */
} Grouplist; /* Type Grouplist is defined */
 

typedef struct hopComm {
    int ngroups;
    int nb;
    float *gdensity;
    float *g1vec;
    float *g2vec;
    float *fdensity;
    Grouplist *gl;
    Slice *s;
} HC;
