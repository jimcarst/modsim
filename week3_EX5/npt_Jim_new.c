#include <stdio.h>
#include <time.h>
#include <math.h>
#include "mt19937.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define NDIM 3
#define N 1000

/* Initialization variables */
const int mc_steps = 100000;
const int output_steps = 1000;
const double packing_fraction = 0.3;
const double diameter = 1.0;
const double delta  = 0.05;
/* Volume change -deltaV, delta V */
const double deltaV = 0.5;
/* Reduced pressure \beta P */
const double betaP = 0.1;
const char* init_filename = "fcc.dat";

/* Simulation variables */
int n_particles = 0;
double radius;
double particle_volume;
double r[N][NDIM];
double box[NDIM];


/* Functions */

/* Return 0 if overlap else 1 if no overlap */
int check_overlap(void){
    int i,j,d;
    int no_overlap = 1;
    for(i = 0; i < n_particles; ++i){
        for(j = i+1; j < n_particles; ++j){
            if(i == j) continue;
            double dist = 0.0;
            for(d = 0; d < NDIM; ++d){
                double min_d = r[j][d] - r[i][d];
                min_d -= (int)(2.0 * min_d / box[d]) * box[d];
                assert(min_d <= 0.5 * box[d]);
                assert(min_d >= -0.5 * box[d]);
                dist += min_d * min_d;
            }
            if(dist <= diameter * diameter){
                no_overlap = 0;
            }
        }
    }
    return no_overlap;
}

int change_volume(void){
    int d,n;
    double volume = 1.0;
    double scale_factor=1.0;
    for(d = 0; d < NDIM; ++d) volume *= box[d];
    assert(volume>0);
    //Change the volume by an amount between (-deltaV,deltaV)
    double newvolume = volume + deltaV*(2.0*dsfmt_genrand() - 1.0);
    assert(newvolume>0);

    //Scale the box and coordinates
    scale_factor = pow(newvolume / volume, 1.0 / NDIM);
    for(n = 0; n < n_particles; ++n){
        for(d = 0; d < NDIM; ++d) r[n][d] *= scale_factor;
    }
    for(d = 0; d < NDIM; ++d) box[d] *= scale_factor;

    //Acceptance criteria for volume change move
    double prob = exp(-betaP*(newvolume-volume)+(n_particles+1)*log(newvolume/volume));
    if(dsfmt_genrand() < prob){
        //Check for overlaps
        if(check_overlap()==1){
            //Accept if no overlaps
            return 1;
        }
    }

    //Reject move and reset to original packing fraction and return
    scale_factor = pow(volume / newvolume, 1.0 / NDIM);
    for(n = 0; n < n_particles; ++n){
        for(d = 0; d < NDIM; ++d) r[n][d] *= scale_factor;
    }
    for(d = 0; d < NDIM; ++d) box[d] *= scale_factor;
    //printf("Reject: %f\n",prob);

    return 0;
}

/* Read the initial configuration provided */
void read_data(void){
    FILE* fp = fopen(init_filename, "r");
    if(fp==NULL){
        printf("Input file not found! Exiting...\n");
        exit(1);
    }
    int n, d, dummy;
    double dmin, dmax, dia;
    dummy = fscanf(fp, "%d\n", &n_particles);
    assert(n_particles<N);
    for(d = 0; d < NDIM; ++d){
        dummy = fscanf(fp, "%lf %lf\n", &dmin, &dmax);
        box[d] = fabs(dmax-dmin);
    }
    for(n = 0; n < n_particles; ++n){
        for(d = 0; d < NDIM; ++d){
           dummy = fscanf(fp, "%lf\t", &r[n][d]);
        }
        dummy = fscanf(fp, "%lf\n", &dia);
    }
    fclose(fp);

    double volume = 1.0;
    for(d = 0; d < NDIM; ++d) volume*= box[d];
    printf("Input: Npart: %d, packing fraction:%lf\n",n_particles,n_particles*M_PI/(6*volume));
}

/* Move a particle randomly */
int move_particle(void){
    //Choose a random particle
    int rpid = n_particles * dsfmt_genrand();
    
    double new_pos[NDIM];
    int d;
    for(d = 0; d < NDIM; ++d){
        //Displace by a random value between -delta and delta
        new_pos[d] = r[rpid][d] + delta * (2.0 * dsfmt_genrand() - 1.0) + box[d];
        //Apply periodic boundaries
        new_pos[d] -= (int)(new_pos[d] / box[d]) * box[d];
        assert(new_pos[d] < box[d]);
        assert(new_pos[d] >= 0.0);
    }

    //Check for overlaps
    int n;
    for(n = 0; n < n_particles; ++n){
        if(n == rpid) continue;
        double dist = 0.0;
        for(d = 0; d < NDIM; ++d){
            double min_d = new_pos[d] - r[n][d];
            // Find the distance with the Nearest Image Convension
            min_d -= (int)(2.0 * min_d / box[d]) * box[d];
            assert(min_d <= 0.5 * box[d]);
            assert(min_d >= -0.5 * box[d]);
            dist += min_d * min_d;
        }
        if(dist <= diameter * diameter){
            //reject the move
            return 0;
        }
    }

    //Accept the move if reached here
    for(d = 0; d < NDIM; ++d) r[rpid][d] = new_pos[d];

    return 1;
}

/* Write the configuration files */
void write_data(int step){
    char buffer[128];
    sprintf(buffer, "coords_step%07d.output", step);
    FILE* fp = fopen(buffer, "w");
    int d, n;
    fprintf(fp, "%d\n", n_particles);
    for(d = 0; d < NDIM; ++d){
        fprintf(fp, "%lf %lf\n",0.0,box[d]);
    }
    for(n = 0; n < n_particles; ++n){
        for(d = 0; d < NDIM; ++d) fprintf(fp, "%lf\t", r[n][d]);
        fprintf(fp, "%lf\n", diameter);
    }
    fclose(fp);
}

/* Scales the box volume and coordinates read from 'init_fileaname'
 * according to the value of 'packing fraction'
 * */
void set_packing_fraction(void){
    double volume = 1.0;
    int d, n;
    for(d = 0; d < NDIM; ++d) volume *= box[d];

    double target_volume = (n_particles * particle_volume) / packing_fraction;
    double scale_factor = pow(target_volume / volume, 1.0 / NDIM);

    for(n = 0; n < n_particles; ++n){
        for(d = 0; d < NDIM; ++d) r[n][d] *= scale_factor;
    }
    for(d = 0; d < NDIM; ++d) box[d] *= scale_factor;
    
    printf("Packing fraction set to: %lf\n",n_particles*M_PI/(6*target_volume));
}

int main(int argc, char* argv[]){

    radius = 0.5 * diameter;

    if(NDIM == 3) particle_volume = M_PI * pow(diameter, 3.0) / 6.0;
    else if(NDIM == 2) particle_volume = M_PI * pow(radius, 2.0);
    else{
        printf("Number of dimensions NDIM = %d, not supported.", NDIM);
        return 0;
    }

    //Read the input configuration
    read_data();

    if(n_particles <= 0){
        printf("Error: Number of particles, n_particles <= 0.\n");
        return 0;
    }

    //Set the packing fraction according to the set variable
    set_packing_fraction();

    dsfmt_seed(time(NULL));
            
    printf("#Step \t Volume \t Move-acceptance\t Volume-acceptance \n");

    int move_accepted = 0;
    int vol_accepted = 0;
    int step, n;
    //Perform MC moves
    for(step = 0; step < mc_steps; ++step){
        for(n = 0; n < n_particles; ++n){
            move_accepted += move_particle();
        }
        //Perform a volume change move
        vol_accepted += change_volume();

        if(step % output_steps == 0){
            printf("%d \t %lf \t %lf \t %lf \n", 
                    step, box[0] * box[1] * box[2], 
                    (double)move_accepted / (n_particles * output_steps), 
                    (double)vol_accepted /  output_steps);
            move_accepted = 0;
            vol_accepted = 0;
            write_data(step);
        }
    }

    return 0;
}
