/* SLICE.C, Daniel Eisenstein, 1997 */
/* Based on a paper by Daniel Eisenstein & Piet Hut,
"HOP: A New Group-Finding Algorithm for N-body Simulations."
See the included documentation or view it at
http://www.sns.ias.edu/~eisenste/hop/hop_doc.html */
 
/* Version 1.0 (12/15/97) -- Original Release */
 
#include "slice.h"
//#include "macros_and_parameters.h"
 
/* ===================== Creating and Destroying Slices ================ */
 
Slice *newslice()
/* Returns a pointer to a newly allocated Slice variable, with internal
variables initialized */
{
    Slice *s;
    s = (Slice *)malloc(sizeof(Slice));
    s->pid = NULL; s->offset = 0;
    s->px = s->py = s->pz = s->vx = s->vy = s->vz = NULL;
    s->ntag = NULL;
    //s->ID = NULL; /* S Skory */
    s->numpart = s->numlist = 0;
    return s;
}
 
void free_tags(Slice *s)
/* Free the tag vector */
{
    if (s->ntag!=NULL) {
      free_ivector(s->ntag, 1, s->numlist);
      s->ntag=NULL;
      //free_ivector(s->ID, 1, s->numlist); /* S Skory */
     //s->ID=NULL;
    }
    return;
}
 
void free_data(Slice *s)
/* Free all the data vectors */
{
    if (s->pid!=NULL) {free(s->pid); s->pid=NULL;}
    if (s->px!=NULL) {free_vector(s->px,1,s->numlist); s->px=NULL;}
    if (s->py!=NULL) {free_vector(s->py,1,s->numlist); s->py=NULL;}
    if (s->pz!=NULL) {free_vector(s->pz,1,s->numlist); s->pz=NULL;}
    if (s->vx!=NULL) {free_vector(s->vx,1,s->numlist); s->vx=NULL;}
    if (s->vy!=NULL) {free_vector(s->vy,1,s->numlist); s->vy=NULL;}
    if (s->vz!=NULL) {free_vector(s->vz,1,s->numlist); s->vz=NULL;}
    return;
}
 
void free_slice(Slice *s)
/* Free the space associated with the vectors in the given Slice */
/* Then free the Slice variable itself */
{
    free_tags(s);
    free_data(s);
    free(s);
    return;
}
 
/* =================================================================== */
 
int f77write(FILE *f, void *p, int len)
/* len is number of bytes to be written from p[0..len-1] */
/* Return 0 if successful, 1 if not */
{
    if (fwrite(&len,sizeof(int),1,f)!=1) return 1;
    if (fwrite(p,1,len,f)!=len) return 1;
    if (fwrite(&len,sizeof(int),1,f)!=1) return 1;
    return 0;
}
 
int f77read(FILE *f, void *p, int maxbytes)
/* Read a FORTRAN style block from the given file */
/* maxbytes is the amount of space the pointer p points to */
/* Space must be allocated to read the whole block into p */
/* Return amount read, scream if there's a problem */
/* Reading is done ZERO-OFFSET */
{
    int size, size2;
    if (fread(&size,sizeof(int),1,f)!=1)
	myerror("f77read(): Error reading begin delimiter.");
    if (size>maxbytes)
	myerror("f77read(): Block delimiter exceeds size of storage.");
    if (size<maxbytes)
	mywarn("f77read(): Block size is smaller than size of storage.");
    if (fread(p,1,size,f)!=size) myerror("f77read(): Error reading data.");
    if (fread(&size2,sizeof(int),1,f)!=1)
	myerror("f77read(): Error reading end delimiter.");
    if (size!=size2) myerror("f77read(): Delimiters do not match.");
    return size;
}
 
/* =================================================================== */
/* The following are public-domain routines from Numerical Repices in C,
2nd edition, by Press, Teulkolsky, Vetterling, & Flannery, 1992, Cambridge
Univ. Press */
 
#define NR_END 1
#define FREE_ARG char*
 
float *vector(long nl, long nh)
/* allocate a float vector with subscript range v[nl..nh] */
{
        float *v;
 
        v=(float *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float)));
        if (!v) myerror("allocation failure in vector()");
        return v-nl+NR_END;
}
 
int *ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
        int *v;
 
        v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
        if (!v) myerror("allocation failure in ivector()");
        return v-nl+NR_END;
}
void free_vector(float *v, long nl, long nh)
/* free a float vector allocated with vector() */
{
        free((FREE_ARG) (v+nl-NR_END));
}
 
void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
        free((FREE_ARG) (v+nl-NR_END));
}
 
#ifdef NOT_USED
/* =================================================================== */
/* =================================================================== */
 
/* Everything from here on down may be expendable */
 
int read_header(FILE *f, Slice *s)
/* Read the header from file f.  Place information and calculated
numbers into the Slice variable, unless *s is NULL, in which case just
skip the header */
/* Return 0 if successful, non-zero otherwise */
{
    int sizeheader[2];
    float header[100];
 
    f77read(f, header, 400);
    f77read(f, sizeheader, 8);
 
    if (s==NULL) return 0;	/* We've been instructed just to skip the
			header without recording the information */
    /* Now read the interesting information out of these arrays */
    s->numpart = sizeheader[0];
    s->numblocks = sizeheader[1];
    s->numperblock = sizeheader[0]/sizeheader[1];
    if (s->numpart != s->numblocks*(s->numperblock))
	myerror("Number of blocks not an even divisor of number of particles.");
 
    s->z = header[0];
    s->boxsize = header[1]*1000.0;	/* We use kpc, not Mpc */
    s->physsize = s->boxsize/(1.0+s->z);	/* We use kpc, not Mpc */
    s->velscale = 100.0*header[1]*sqrt(3.0/8.0/PI)/(1.0+s->z);
		/* To go from data to pec vel */
    s->omega = header[4];
    if (header[6]!=0.0) myerror("HDM component listed in header.");
    s->lambda = header[7];
    s->h0 = header[8];
    s->sigma8 = header[9];	/* At z=0 */
 
    /* Now find some computed quantities. */
    s->a = 1.0/(1.0+s->z);
    s->curv = 1.0-s->omega-s->lambda;
    s->gamma = s->omega*(s->h0);
    s->specn = 1.0;
    s->hubb = 0.1*sqrt(s->omega/CUBE(s->a)+s->curv/SQR(s->a)+s->lambda)*(s->a);
 
    /* The following assume Omega = 1 */
    s->masspart = RHOCRIT/s->numpart*CUBE(s->boxsize);
    s->growth = s->a;
    s->t = HTIME*(s->h0)*pow(s->a, 1.5);
    return 0;
}
 
void normalizedata(Slice *s, int conp, int conv)
/* Put raw data into comoving h^-1 kpc and km/s units */
{
    int j;
    float velnorm;
    if (conp) {
	for (j=1;j<=s->numlist;j++) s->px[j] *= s->boxsize;
	for (j=1;j<=s->numlist;j++) s->py[j] *= s->boxsize;
	for (j=1;j<=s->numlist;j++) s->pz[j] *= s->boxsize;
    }
 
    if (conv) {
	for (j=1;j<=s->numlist;j++) s->vx[j] *= s->velscale;
	for (j=1;j<=s->numlist;j++) s->vy[j] *= s->velscale;
	for (j=1;j<=s->numlist;j++) s->vz[j] *= s->velscale;
    }
    return;
}
 
/* ================================================================ */
 
int read_alldata(FILE *f, FILE *ftag, Slice *s, int conp, int conv)
/* Read all the data, including the tags if ftag!=NULL. */
/* Store positions and velocities unless conp or conv = 0 */
/* Assume that the data file header has been skipped, but read the
tag file header. */
{
    int block;
    float *dummylist;
 
    dummylist = NULL;
    if (!conp || !conv) dummylist = vector(1,s->numperblock);
    if (s->pid!=NULL)
	mywarn("Non-NULL s->pid[] passed to read_alldata().  Ignoring...");
 
    s->numlist=s->numpart;
    if (conp) {
	s->px=vector(1,s->numlist);
	s->py=vector(1,s->numlist);
	s->pz=vector(1,s->numlist);
    }
    if (conv) {
	s->vx=vector(1,s->numlist);
	s->vy=vector(1,s->numlist);
	s->vz=vector(1,s->numlist);
    }
    if (ftag!=NULL) {
      s->ntag = ivector(1,s->numlist);
       s->ID = ivector(1,s->numlist); /* S Skory */
    }
    
 
    printf("Reading data...");
    for (block=0;block<s->numblocks;block++) {
	/* Read the three position blocks */
	if (conp) {  /* Store the data */
	    f77read(f, s->px+s->numperblock*block+1, s->numperblock*sizeof(float));
	    f77read(f, s->py+s->numperblock*block+1, s->numperblock*sizeof(float));
	    f77read(f, s->pz+s->numperblock*block+1, s->numperblock*sizeof(float));
	} else {     /* Don't store the data */
	    f77read(f, dummylist+1, s->numperblock*sizeof(float));
	    f77read(f, dummylist+1, s->numperblock*sizeof(float));
	    f77read(f, dummylist+1, s->numperblock*sizeof(float));
	}
	/* Now read the three velocity blocks */
	if (conv) {  /* Store the data */
	    f77read(f, s->vx+s->numperblock*block+1, s->numperblock*sizeof(float));
	    f77read(f, s->vy+s->numperblock*block+1, s->numperblock*sizeof(float));
	    f77read(f, s->vz+s->numperblock*block+1, s->numperblock*sizeof(float));
	} else {     /* Don't store the data */
	    f77read(f, dummylist+1, s->numperblock*sizeof(float));
	    f77read(f, dummylist+1, s->numperblock*sizeof(float));
	    f77read(f, dummylist+1, s->numperblock*sizeof(float));
	}
	if (block%8==1) {printf("."); fflush(stdout);}
    }
    if (dummylist!=NULL) free_vector(dummylist, 1, s->numperblock);
    normalizedata(s,conp,conv);
 
    if (ftag!=NULL) {
	printf("tags..."); fflush(stdout);
	readalltags(ftag, s);
    }
 
    printf("done!"); fflush(stdout);
    return 0;
}
 
int read_partdata(FILE *f, FILE *ftag, Slice *s)
/* Read one block (128k particles) of data into Slice s.  Allocate needed
storage, erasing and freeing previous storage. */
/* This cannot be done with s->pid!=NULL, so s->pid is ignored and s->numlist
is reset to BLOCKSIZE */
/* Unlike other routines, this stores both positions and velocities in
all cases (since the storage requirements are already small */
/* If ftag==NULL, don't read the tag file.  Otherwise do read it. */
{
    if (s->pid!=NULL)
	mywarn("Non-trivial pid[] not supported with incremental reads");
    /* If we need to reallocate memory, do it.  Otherwise, just write over */
    if (s->px==NULL || s->vx==NULL || s->numlist!=s->numperblock) {
	if (s->px!=NULL) free_vector(s->px,1,s->numlist);
	if (s->py!=NULL) free_vector(s->py,1,s->numlist);
	if (s->pz!=NULL) free_vector(s->pz,1,s->numlist);
	if (s->vx!=NULL) free_vector(s->vx,1,s->numlist);
	if (s->vy!=NULL) free_vector(s->vy,1,s->numlist);
	if (s->vz!=NULL) free_vector(s->vz,1,s->numlist);
	if (ftag!=NULL && s->ntag!=NULL) {
	  free_ivector(s->ntag, 1, s->numlist);
	  free_ivector(s->ID, 1, s->numlist); /* S Skory */
	}
	s->numlist = s->numperblock;
	s->px = vector(1,s->numlist);
	s->py = vector(1,s->numlist);
	s->pz = vector(1,s->numlist);
	s->vx = vector(1,s->numlist);
	s->vy = vector(1,s->numlist);
	s->vz = vector(1,s->numlist);
	if (ftag!=NULL) {
	  s->ntag = ivector(1, s->numlist);
	  s->ID = ivector(1, s->numlist); /* S Skory */
	}
	s->offset=0;
	/* fprintf(stderr, "Reallocating data arrays.\n"); */
    }
    else s->offset+=s->numlist;
 
    f77read(f,s->px+1,sizeof(float)*s->numlist);
    f77read(f,s->py+1,sizeof(float)*s->numlist);
    f77read(f,s->pz+1,sizeof(float)*s->numlist);
    f77read(f,s->vx+1,sizeof(float)*s->numlist);
    f77read(f,s->vy+1,sizeof(float)*s->numlist);
    f77read(f,s->vz+1,sizeof(float)*s->numlist);
 
    if (ftag!=NULL) readtag(ftag, s->numlist, s->ntag);
    normalizedata(s,1,1);
    return 0;
}
 
/* =============================================================== */
 
int readtag(FILE *f, int numread, int *ntag)
/* Read numread values from FILE f and put the values in ntag */
/* Return 0 if successful, 1 if not */
/* The storage ntag[1..numread] must exist */
/* Note: the first 8 bytes of the tag file contain the number of particles
and the number of groups.  These must be skipped before calling this routine. */
{
    if (fread(ntag+1, sizeof(int), numread, f)!=numread)
	myerror("Error in reading tag file.");
    return 0;
}
 
int skiptagheader(FILE *f, Slice *s)
/* Read the first 8 bytes from the tag file.  Check that the first int equals
the number of particles.  Return the second, which is the number of groups */
{
    int dummy[2];
    if (fread(&dummy, sizeof(int), 2, f)!=2) myerror("Error in reading tag file.");
    if (s->numpart!=0 && dummy[0]!=s->numpart)
	myerror("First number in tag file doesn't match expected number of particles.");
    s->numgroups = dummy[1];
    return dummy[1];
}
 
int readalltags(FILE *f, Slice *s)
/* Read the whole tag file.  Allocate memory as needed */
/* Return the number of groups */
{
    int dummy[2];
    if (s->ntag==NULL || s->numlist!=s->numpart) {
	if (s->ntag!=NULL) {
	  free_ivector(s->ntag, 1, s->numlist);
	  free_ivector(s->ID, 1, s->numlist); /* S Skory */
	}
	s->numlist = s->numpart;
	s->ntag = ivector(1, s->numlist);
	s->ID = ivector(1, s->numlist);
    }
    if (fread(&dummy, sizeof(int), 2, f)!=2) myerror("Error 1 in reading tag file.");
    if (dummy[0]!=s->numpart)
	myerror("First int of tag file doesn't match numpart.");
    s->numgroups = dummy[1];
 
    if (fread(s->ntag+1, sizeof(int), s->numlist, f)!=s->numlist)
	myerror("Couldn't read entire tag file.");
    return dummy[1];
}
 
#endif
 
/* ===================== Warnings and Errors =========================== */
 
/* Print a message and die */
void myerror(char *message)
{
    fprintf(stderr, "%s\n", message);
    exit(1); return;
}
 
/* Just print a message */
void mywarn(char *message)
{
    fprintf(stderr, "%s\n", message);
    fflush(NULL);  /* Flush everything, so we know where we are */
    return;
}
