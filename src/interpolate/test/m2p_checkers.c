#include "m2p_checkers.h"

static char stamp_checker[] = "***\nmodule '' __FILE__ ''\ncompiled '' __TIMESTAMP__ ''\n***";

void createRandomParticlePosition_2d(
			        int   p_per_cell, 
			        float *physmin, 
	       			float *physmax,
					int   *mesh_size,
		    		float *xp,
					int   ndim
){
	int k, x, y, z, p;
	float h[2];

	srand(0);

	h[0] = (physmax[0] - physmin[0])/(mesh_size[0] - 1);
	h[1] = (physmax[1] - physmin[1])/(mesh_size[1] - 1);
	
	p = 0;
	for(k = 0; k < p_per_cell; k++){
		for(y = 2; y < mesh_size[1] - 2; y++){
			for(x = 2; x < mesh_size[0] - 2; x++){
				xp[ndim*p]     = (x*h[0] + physmin[0]);
				xp[ndim*p + 1] = (y*h[1] + physmin[1]);

				xp[ndim*p]   += ((rand()%10000 - 5000.0f)/10000.0f)*h[0];
				xp[ndim*p+1] += ((rand()%10000 - 5000.0f)/10000.0f)*h[1];
				p++;
			}
		}
	}
}

void createRandomParticlePosition_3d(
			        int   p_per_cell, 
			        float *physmin, 
	       			float *physmax,
					int   *mesh_size,
		    		float *xp,
					int    ndim
){
	int k, x, y, z, p;
	float h[3];

	srand(0);

	h[0] = (physmax[0] - physmin[0])/(mesh_size[0] - 1);
	h[1] = (physmax[1] - physmin[1])/(mesh_size[1] - 1);
	h[2] = (physmax[2] - physmin[2])/(mesh_size[2] - 1);
	
	p = 0;
	for(k = 0; k < p_per_cell; k++){
		for(z = 2; z < mesh_size[2] - 2; z++){
			for(y = 2; y < mesh_size[1] - 2; y++){
				for(x = 2; x < mesh_size[0] - 2; x++){
					xp[ndim*p]     = (x*h[0] + physmin[0]);
					xp[ndim*p + 1] = (y*h[1] + physmin[1]);
					xp[ndim*p + 2] = (z*h[2] + physmin[2]);

					xp[ndim*p]   += ((rand()%10000 - 5000.0f)/10000.0f)*h[0];
					xp[ndim*p+1] += ((rand()%10000 - 5000.0f)/10000.0f)*h[1];
					xp[ndim*p+2] += ((rand()%10000 - 5000.0f)/10000.0f)*h[2];
					p++;
				}
			}
		}
	}
}

void createRandomParticleMass_2d(
			        int   p_per_cell, 
					int   *mesh_size,
                    float *mass
){
	int k, x, y, z, p;
	float h[3];

	srand(0);

	p = 0;
	for(k = 0; k < p_per_cell; k++){
		for(y = 0; y < mesh_size[1]; y++){
			for(y = 0; y < mesh_size[1]; y++){
				for(x = 0; x < mesh_size[0]; x++){
					mass[p] = ((rand()%10001 - 5000)/10000.0);
					p++;
				}
			}
		}
	}
}


void createRandomParticleMass_3d(
			        int   p_per_cell, 
					int   *mesh_size,
                    float *mass
){
	int k, x, y, z, p;
	float h[3];

	srand(0);

	p = 0;
	for(k = 0; k < p_per_cell; k++){
		for(z = 0; z < mesh_size[2]; z++){
			for(y = 0; y < mesh_size[1]; y++){
				for(y = 0; y < mesh_size[1]; y++){
					for(x = 0; x < mesh_size[0]; x++){
						mass[p] = ((rand()%10001 - 5000)/10000.0);
						p++;
					}
				}
			}
		}
	}
}

void createRandomMesh_2d_m2p(
					int   *mesh_size,
					int    nmass,
					float *mesh
){
	int i, x, y, z, mesh_ID;

	srand(0);

	for(y = 0; y < mesh_size[1]; y++){
		for(x = 0; x < mesh_size[0]; x++){
			mesh_ID = x + y*mesh_size[0];

			for(i = 0; i < nmass; i++){
				mesh[nmass*mesh_ID + i] = ((rand()%10001 - 5000)/10000.0);
			}
		}
	}
}

void createRandomMesh_3d_m2p(
					int   *mesh_size,
					int    nmass,
					float *mesh
){
	int i, x, y, z, mesh_ID;

	srand(0);

	for(z = 0; z < mesh_size[2]; z++){
		for(y = 0; y < mesh_size[1]; y++){
			for(x = 0; x < mesh_size[0]; x++){
				mesh_ID = x + y*mesh_size[0] + z*mesh_size[0]*mesh_size[1];

				for(i = 0; i < nmass; i++){
					mesh[nmass*mesh_ID + i] = ((rand()%10001 - 5000)/10000.0);
				}
			}
		}
	}
}

void check_interpolation_2d_m2p_mp4(
                             float* xp,
                             float* mass,
                             float* grid,
                             float* physmin,
                             float* physmax,
                             int*   mesh_size,
                             int*   np,
							 int*   ndim_in,
							 int*   nmass_in
) {
    int *offset;
    float *xp_sorted;
    float* grid_seq;
    float* mass_seq;
    float  h[3];
    int    ngrid;
    int   p_id, g_id;
    int x, y, z, p;
    int m_x, m_y, c_x, c_y, cell;
    float m_coor[2], c_x_coor, c_y_coor, c_z_coor;
    float pos_x, pos_y, p_mass;
    cl_int error;
    int ncell;
    int dim;
    int g_id2;
    int *p_cell;
	int  c_ID, m_ID, m[2];
	int  dist[2];

	int mesh_mass;

    int nparticle = *np;
	int ndim 	  = *ndim_in;
	int nmass     = *nmass_in;

    int n = 0;
    float sum = 0;

    float x01, x02, x03, xp1, xp2, xp3, wx1, wx2, wx3, x1, x2, x3;
    int ip1, ip2, ip3;

    ncell = (cell_size_ext[0]);
    ngrid = mesh_size[0];
    for(dim = 1; dim < ndim; dim++){
        ncell *= (cell_size_ext[dim]);
        ngrid *= mesh_size[dim];
    }

    for(dim = 0; dim < ndim; dim++)
        h[dim] = (physmax[dim] - physmin[dim])/(mesh_size[dim] - 1);

    mass_seq = (float *)malloc(sizeof(float) * nparticle * nmass);

    for(p_id = 0; p_id < nparticle*nmass; p_id++)
        mass_seq[p_id] = 0.0f;
    
    for(p = 0; p < nparticle; p++){
        x01 = (xp[p*ndim] - physmin[0])/h[0];
        x02 = (xp[p*ndim + 1] - physmin[1])/h[1];
    
	    ip1 = (int)x01 ;
	    ip2 = (int)x02 ;

	    xp1 = x01 - (int)x01;
	    xp2 = x02 - (int)x02;

        for (y = -1; y < 3; y++){
	    	x2 = fabs(xp2 - (float)y);

			wx2 = 0.0f;
	        if(x2 < 1.0f)
		        wx2 = 1.0f - x2*x2*(2.5f - 1.5f*x2);
	        else if(x2 < 2.0f)
	            wx2= 2.0f + (-4.0f + (2.5f - 0.5f*x2)*x2)*x2;

            for (x = -1; x < 3; x++){
		       x1 = fabs(xp1 - (float)x);

				wx1 = 0.0f;
	    	    if(x1 < 1.0f)
		            wx1 = 1.0f - x1*x1*(2.5f - 1.5f*x1);
	    	    else if(x1 < 2.0f) 
	                wx1 = 2.0f + (-4.0f + (2.5f - 0.5f*x1)*x1)*x1;

                c_x_coor = (ip1 + x)*h[0] + physmin[0];  
                c_y_coor = (ip2 + y)*h[1] + physmin[1]; 

                cell     = (ip1 + x) + (ip2 + y)*mesh_size_ext[0];

				for(mesh_mass = 0; mesh_mass < nmass; mesh_mass++){
                	mass_seq[nmass*p + mesh_mass] += grid[nmass*cell + mesh_mass]*wx1*wx2;
				}
			}

		}
    }

    for(p_id = 0; p_id < nparticle*nmass; p_id++){
        if(fabs(mass[p_id] - mass_seq[p_id]) > 0.01){ // HERE IS IN QUESTION!
            printf("Interpolation is wrong at p_id: %d \n", (p_id/nmass));
            printf("mass[%d]: %12.8f, mass_seq[%d]: %12.8f \n", p_id, mass[p_id], p_id, mass_seq[p_id]);

            for(g_id2 = 0; g_id2 <= 10 + p_id; g_id2++){
                printf("grid[%d]: %f grid_seq[%d]: %f \n", g_id2, mass[g_id2], g_id2, mass_seq[g_id2]);
            }
            exit(0);
        }
    }

    printf("Interpolation results are correct! \n");

    return;
}
void check_interpolation_3d_m2p_mp4(
                             float* xp,
                             float* mass,
                             float* grid,
                             float* physmin,
                             float* physmax,
                             int*   mesh_size,
                             int*   np,
							 int*   ndim_in,
							 int*   nmass_in
) {
    int *offset;
    float *xp_sorted;
    float* grid_seq;
    float* mass_seq;
    float  h[3];
    int    ngrid;
    int   p_id, g_id;
    int x, y, z, p;
    int m_x, m_y, c_x, c_y, cell;
    float m_coor[2], c_x_coor, c_y_coor, c_z_coor;
    float pos_x, pos_y, p_mass;
    cl_int error;
    int ncell;
    int dim;
    int g_id2;
    int *p_cell;
	int  c_ID, m_ID, m[2];
	int  dist[2];

	int mesh_mass;
	
    int nparticle = *np;
	int ndim 	  = *ndim_in;	
	int nmass 	  = *nmass_in;	

    int n = 0;
    float sum = 0;

    float x01, x02, x03, xp1, xp2, xp3, wx1, wx2, wx3, x1, x2, x3;
    int ip1, ip2, ip3;

    ncell = (cell_size_ext[0]);
    ngrid = mesh_size[0];
    for(dim = 1; dim < ndim; dim++){
        ncell *= (cell_size_ext[dim]);
        ngrid *= mesh_size[dim];
    }

    for(dim = 0; dim < ndim; dim++)
        h[dim] = (physmax[dim] - physmin[dim])/(mesh_size[dim] - 1);

    mass_seq = (float *)malloc(sizeof(float) * nparticle * nmass);

    for(p_id = 0; p_id < nparticle*nmass; p_id++)
        mass_seq[p_id] = 0.0f;
    
    for(p = 0; p < nparticle; p++){
        x01 = (xp[p*ndim] - physmin[0])/h[0];
        x02 = (xp[p*ndim + 1] - physmin[1])/h[1];
        x03 = (xp[p*ndim + 2] - physmin[2])/h[2];
    
	    ip1 = (int)x01 ;
	    ip2 = (int)x02 ;
	    ip3 = (int)x03 ;

	    xp1 = x01 - (int)x01;
	    xp2 = x02 - (int)x02;
	    xp3 = x03 - (int)x03;

        for (z = -1; z < 3; z++){
	        x3 = fabs(xp3 - (float)z);

			wx3 = 0.0f;
	        if(x3 < 1.0f)
	            wx3 = 1.0f - x3*x3*(2.5f - 1.5f*x3);
	        else if(x3 < 2.0f)
	            wx3= 2.0f + (-4.0f + (2.5f - 0.5f*x3)*x3)*x3;

            for (y = -1; y < 3; y++){
	            x2 = fabs(xp2 - (float)y);

				wx2 = 0.0f;
	            if(x2 < 1.0f)
		            wx2 = 1.0f - x2*x2*(2.5f - 1.5f*x2);
	            else if(x2 < 2.0f)
	                wx2= 2.0f + (-4.0f + (2.5f - 0.5f*x2)*x2)*x2;

                for (x = -1; x < 3; x++){
		            x1 = fabs(xp1 - (float)x);

					wx1 = 0.0f;
	    	        if(x1 < 1.0f)
		                wx1 = 1.0f - x1*x1*(2.5f - 1.5f*x1);
	    	        else if(x1 < 2.0f) 
	                    wx1 = 2.0f + (-4.0f + (2.5f - 0.5f*x1)*x1)*x1;

                    c_x_coor = (ip1 + x)*h[0] + physmin[0];  
                    c_y_coor = (ip2 + y)*h[1] + physmin[1]; 
                    c_z_coor = (ip3 + z)*h[2] + physmin[2]; 

                    cell     = (ip1 + x) + (ip2 + y)*mesh_size_ext[0] + (ip3 + z)*mesh_size_ext[0]*mesh_size_ext[1];

	
					for(mesh_mass = 0; mesh_mass < nmass; mesh_mass++){
					if(p == 0){
//						printf("%f \n", grid[nmass*cell + mesh_mass]);
					}
   		             	mass_seq[nmass*p + mesh_mass] += grid[nmass*cell + mesh_mass]*wx1*wx2*wx3;
					}
				}
            }
        }        
    }

    for(p_id = 0; p_id < nparticle*nmass; p_id++){
        if(fabs(mass[p_id] - mass_seq[p_id]) > 0.01){ // HERE IS IN QUESTION!
            printf("Interpolation is wrong at p_id: %d \n", (int)(p_id));
            printf("mass[p_id]: %12.8f, mass_seq[p_id]: %12.8f \n", mass[p_id], mass_seq[p_id]);

            for(g_id2 = 0; g_id2 <= p_id+10; g_id2++){
                printf("grid[%d]: %f grid_seq[%d]: %f \n", g_id2, mass[g_id2], g_id2, mass_seq[g_id2]);
            }
            exit(0);
        }
    }

    printf("Interpolation results are correct! \n");

    return;
}

void check_interpolation_2d_m2p_bsp2(
                             float* xp,
                             float* mass,
                             float* grid,
                             float* physmin,
                             float* physmax,
                             int*   mesh_size,
                             int*   np,
							 int*   ndim_in,
							 int*   nmass_in
) {
    int *offset;
    float *xp_sorted;
    float* grid_seq;
    float* mass_seq;
    float  h[3];
    int    ngrid;
    int   p_id, g_id;
    int x, y, z, p;
    int m_x, m_y, c_x, c_y, cell;
    float m_coor[2], c_x_coor, c_y_coor, c_z_coor;
    float pos_x, pos_y, p_mass;
    cl_int error;
    int ncell;
    int dim;
    int g_id2;
    int *p_cell;
	int  c_ID, m_ID, m[2];
	int  dist[2];

	int mesh_mass;

    int nparticle = *np;
    int ndim      = *ndim_in;
    int nmass     = *nmass_in;

    int n = 0;
    float sum = 0;

    float x01, x02, x03, xp1, xp2, xp3, wx1, wx2, wx3, x1, x2, x3;
    int ip1, ip2, ip3;

    ncell = (cell_size_ext[0]);
    ngrid = mesh_size[0];
    for(dim = 1; dim < ndim; dim++){
        ncell *= (cell_size_ext[dim]);
        ngrid *= mesh_size[dim];
    }

    for(dim = 0; dim < ndim; dim++)
        h[dim] = (physmax[dim] - physmin[dim])/(mesh_size[dim] - 1);

    mass_seq = (float *)malloc(sizeof(float) * nparticle * nmass);

    for(p_id = 0; p_id < nparticle*nmass; p_id++)
        mass_seq[p_id] = 0.0f;
    
    for(p = 0; p < nparticle; p++){
        x01 = (xp[p*ndim] - physmin[0])/h[0];
        x02 = (xp[p*ndim + 1] - physmin[1])/h[1];
    
	    ip1 = (int)x01; //x-origin of the cell
	    ip2 = (int)x02; //y-origin of the cell

	    xp1 = x01 - (int)x01; // h-distance in x
	    xp2 = x02 - (int)x02; // h-distance in y

        for (y = 0; y < 2; y++){
	    	x2 = fabs(xp2 - (float)y);
			wx2 = 1.0f - x2;

            for (x = 0; x < 2; x++){
		        x1 = fabs(xp1 - (float)x);
				wx1 = 1.0f - x1;

                cell     = (ip1 + x) + (ip2 + y)*mesh_size_ext[0];

				for(mesh_mass = 0; mesh_mass < nmass; mesh_mass++){
                	mass_seq[nmass*p + mesh_mass] += grid[nmass*cell + mesh_mass]*wx1*wx2;
				}
			}

		}
    }

    for(p_id = 0; p_id < nparticle*nmass; p_id++){
        if(fabs(mass[p_id] - mass_seq[p_id]) > 0.01){ // HERE IS IN QUESTION!
            printf("Interpolation is wrong at p_id: %d \n", (p_id/nmass));
            printf("mass[%d]: %12.8f, mass_seq[%d]: %12.8f \n", p_id, mass[p_id], p_id, mass_seq[p_id]);

            for(g_id2 = 0; g_id2 <= 10 + p_id; g_id2++){
                printf("grid[%d]: %f grid_seq[%d]: %f \n", g_id2, mass[g_id2], g_id2, mass_seq[g_id2]);
            }
            exit(0);
        }
    }

    printf("Interpolation results are correct! \n");

    return;
}

void check_interpolation_3d_m2p_bsp2(
                             float* xp,
                             float* mass,
                             float* grid,
                             float* physmin,
                             float* physmax,
                             int*   mesh_size,
                             int*   np,
							 int*   ndim_in,
							 int*   nmass_in
) {
    int *offset;
    float *xp_sorted;
    float* grid_seq;
    float* mass_seq;
    float  h[3];
    int    ngrid;
    int   p_id, g_id;
    int x, y, z, p;
    int m_x, m_y, c_x, c_y, cell;
    float m_coor[2], c_x_coor, c_y_coor, c_z_coor;
    float pos_x, pos_y, p_mass;
    cl_int error;
    int ncell;
    int dim;
    int g_id2;
    int *p_cell;
	int  c_ID, m_ID, m[2];
	int  dist[2];

	int mesh_mass;

    int nparticle = *np;
	int ndim      = *ndim_in;
	int nmass     = *nmass_in;

    int n = 0;
    float sum = 0;

    float x01, x02, x03, xp1, xp2, xp3, wx1, wx2, wx3, x1, x2, x3;
    int ip1, ip2, ip3;

    ncell = (cell_size_ext[0]);
    ngrid = mesh_size[0];
    for(dim = 1; dim < ndim; dim++){
        ncell *= (cell_size_ext[dim]);
        ngrid *= mesh_size[dim];
    }

    for(dim = 0; dim < ndim; dim++)
        h[dim] = (physmax[dim] - physmin[dim])/(mesh_size[dim] - 1);

    mass_seq = (float *)malloc(sizeof(float) * nparticle * nmass);

    for(p_id = 0; p_id < nparticle*nmass; p_id++)
        mass_seq[p_id] = 0.0f;
    
    for(p = 0; p < nparticle; p++){
        x01 = (xp[p*ndim] - physmin[0])/h[0];
        x02 = (xp[p*ndim + 1] - physmin[1])/h[1];
        x03 = (xp[p*ndim + 2] - physmin[2])/h[2];
    
	    ip1 = (int)x01; //x-origin of the cell
	    ip2 = (int)x02; //y-origin of the cell
	    ip3 = (int)x03; //z-origin of the cell

	    xp1 = x01 - (int)x01; // h-distance in x
	    xp2 = x02 - (int)x02; // h-distance in y
	    xp3 = x03 - (int)x03; // h-distance in y

        for (z = 0; z < 2; z++){
	    	x3 = fabs(xp3 - (float)z);
			wx3 = 1.0f - x3;

	        for (y = 0; y < 2; y++){
		    	x2 = fabs(xp2 - (float)y);
				wx2 = 1.0f - x2;

          	    for (x = 0; x < 2; x++){
		    	    x1 = fabs(xp1 - (float)x);
					wx1 = 1.0f - x1;

                	cell     = (ip1 + x) + (ip2 + y)*mesh_size_ext[0] + (ip3 + z)*mesh_size_ext[0]*mesh_size_ext[1];

					for(mesh_mass = 0; mesh_mass < nmass; mesh_mass++){
                		mass_seq[nmass*p + mesh_mass] += grid[nmass*cell + mesh_mass]*wx1*wx2*wx3;
					}
				}
			}
		}
    }

    for(p_id = 0; p_id < nparticle*nmass; p_id++){
        if(fabs(mass[p_id] - mass_seq[p_id]) > 0.01){ // HERE IS IN QUESTION!
            printf("Interpolation is wrong at p_id: %d \n", (p_id/nmass));
            printf("mass[%d]: %12.8f, mass_seq[%d]: %12.8f \n", p_id, mass[p_id], p_id, mass_seq[p_id]);

            for(g_id2 = 0; g_id2 <= 10 + p_id; g_id2++){
                printf("grid[%d]: %f grid_seq[%d]: %f \n", g_id2, mass[g_id2], g_id2, mass_seq[g_id2]);
            }
            exit(0);
        }
    }

    printf("Interpolation results are correct! \n");

    return;
}

