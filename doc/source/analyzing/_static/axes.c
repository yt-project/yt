#include "axes.h"

void calculate_axes(ParticleCollection *part,
    double *ax1, double *ax2, double *ax3)
{
    int i;
    for (i = 0; i < part->npart; i++) {
        if (ax1[0] > part->xpos[i]) ax1[0] = part->xpos[i];
        if (ax2[0] > part->ypos[i]) ax2[0] = part->ypos[i];
        if (ax3[0] > part->zpos[i]) ax3[0] = part->zpos[i];
        if (ax1[1] < part->xpos[i]) ax1[1] = part->xpos[i];
        if (ax2[1] < part->ypos[i]) ax2[1] = part->ypos[i];
        if (ax3[1] < part->zpos[i]) ax3[1] = part->zpos[i];
    }
}
