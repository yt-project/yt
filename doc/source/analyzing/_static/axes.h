typedef struct structParticleCollection {
     long npart;
     double *xpos;
     double *ypos;
     double *zpos;
} ParticleCollection;

void calculate_axes(ParticleCollection *part,
         double *ax1, double *ax2, double *ax3);
