






// Version 2 of the RNG
//-----------------------------------
// private per-thread generator state.  
unsigned int TausStep(unsigned int *z, int S1, int S2, int S3, unsigned long M)  
{  
	unsigned b=(((*z << S1) ^ *z) >> S2);
	*z = (((*z & M) << S3) ^ b);
	return *z;
}  


// A and C are constants  
unsigned int LCGStep(unsigned int *z, unsigned int A, unsigned int C)  
{  
	*z=(A * (*z)+C);
	return *z;
} 


float HybridTaus(unsigned int *z1, unsigned int *z2, unsigned int *z3, unsigned int *z4)  
{  
	// Combined period is lcm(p1,p2,p3,p4)~ 2^121
	return 2.3283064365387e-10 * (                                          // Periods
			TausStep(z1, 13, 19, 12, 4294967294UL) ^  // p1=2^31-1
			TausStep(z2, 2, 25, 4, 4294967288UL) ^    // p2=2^30-1
			TausStep(z3, 3, 11, 17, 4294967280UL) ^   // p3=2^28-1
			LCGStep(z4, 1664525, 1013904223UL)        // p4=2^32
	);
}  

float rand(unsigned int *z1, unsigned int *z2, unsigned int *z3, unsigned int *z4)
{
	return HybridTaus(z1, z2, z3, z4);
}




//===============================================================
// Tissue structure.
//------------------
typedef struct {
	float mu_a;     // Absorption coefficient
	float mu_s;     // Scattering coefficient
	float g;        // Anisotropy
	float albedo;
} Tissue;





//===============================================================
// The photon structure.
//----------------------
#define ALIVE 1
#define DEAD  0

typedef struct {
	// Location of the photon in carteesian coordinates.
	float x;
	float y;
	float z;

	// Direction of photon.
	float ux;
	float uy;
	float uz;

	// Weight of the photon.
	float weight;

	// Step size of the photon
	float step_size;

	// Boolean tracking whether photon is alive or dead.
	int status;

} Photon;




//==============================================================
// Initialize the photon's attributes before propagation through
// the medium commences.
//-----------------------
void InitPhoton(Photon *p)
{
	// Initialize the starting coordinates.
	p->x = 0.0f;
	p->y = 0.0f;
	p->z = 0.0f;

	// Set the initial direction of the photon.
	p->ux = 0.0f;
	p->uy = 0.0f;
	p->uz = 1.0f;

	// Initialize the starting weight of the photon.
	p->weight = 1.0f;

	// Initialize the step size.
	p->step_size = 0.0f;

	// Photon was just created, therefore we set it to ALIVE.
	p->status = ALIVE;


}



//==============================================================
// Initialize the values for the tissue (i.e. scattering,
// absorption, etc.).
//-------------------
void InitTissue(Tissue *t)
{
	t->mu_a = 1.0f;  // cm^-1
	t->mu_s = 100.0f; // cm^-1
	t->g    = 0.9f;
	t->albedo = t->mu_s / (t->mu_s + t->mu_a);
}





//void testRng(const int NUM_RANDS, unsigned int *z1, unsigned int *z2,
//             unsigned int *z3, unsigned int *z4, __global float *results)
//{
//    int i;
//    int index = get_global_id(0);
//    int g_work_size = get_global_size(0);
//    for (i = 0; i < NUM_RANDS; i++)
//        results[index*NUM_RANDS + i] = rand(&z1, &z2, &z3, &z4);
//
//}

//#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics : enable;

// NOT SUPPORTED ON CURRENT GPU!!!
//#pragma OPENCL EXTENSION cl_amd_fp64 : enable;

//==============================================================
// Main execution kernel.
//-----------------------
__kernel void LaunchPhoton(__global float *initial_state_vals,
		__global float *global_results,
		const int NUMPHOTONS,
		const int DETECTOR_SIZE)
{

	// Get indices and size of the work-items.
	int index = get_global_id(0);
	int global_work_size = get_global_size(0);


	// Zero out the results array since we will be accumulating values.
	int k;
	for (k = 0; k < DETECTOR_SIZE * global_work_size; k++) {
		global_results[k] = 0;
	}


	const float CHANCE = 0.1f;     // Chance of surviving roulette.
	const float THRESHOLD = 0.01f; // Threshold for determining of we should perform roulette.


	


	// Each work-item (i.e. thread) needs it's own RNG.  z1 - z4 are used
	// as state for each work-item's RNG.
	unsigned int z1;
	unsigned int z2;
	unsigned int z3;
	unsigned int z4;

	// Give each work-item it's 4 random seed values, which will also
	// hold the state of the RNG.
	z1 = initial_state_vals[index*4 + 0];
	z2 = initial_state_vals[index*4 + 1];
	z3 = initial_state_vals[index*4 + 2];
	z4 = initial_state_vals[index*4 + 3];


	// Temporary photon trajectory directions.
	float uxx, uyy, uzz;


	// Test the RNG briefly.
	//testRng(NUMPHOTONS, &z1, &z2, &z3, &z4, global_results);



	// FIXME:  These should not be defined as constants.
	//==================================================
	//          Constants defined for our detection geometry.
	float radial_size 	= 3.0;   /* cm, total range over which bins extend */
	int NR          	= 100;	 /* set number of bins.  */
	/* IF NR IS ALTERED, THEN USER MUST ALSO ALTER THE ARRAY DECLARATION TO A SIZE = NR + 1. */
	float dr = radial_size/(float)NR;  /* cm */


	// For SPIN portion.
	float cost, sint;	// cosine and sine of the polar deflection angle theta.
	float cosp, sinp;	// cosine and sine of the azimuthal angle psi.
	float psi;

	// Self explanatory.
	const float PI = 3.14159265;



	//---------------------------------------------------------
	// Create the photon and initialize it.
	Photon p;
	//InitPhoton(&p);


	//---------------------------------------------------------
	// Create the tissue and initialize it.
	Tissue t;
	InitTissue(&t);



	//---------------------------------------------------------
	// Perform Hop, Drop, Spin, and Roulette.
	int i = 0;
	int ir = 0;
	float r = 0.0f;
	float absorbed = 0.0f;

	//
	// FIXME: Should not be hard coded value here.
	//--------------------------------------------
	float g = t.g;
	float tempdir = 0.0f;
	//float uxx, uyy, uzz;




	int j = 0;
	const float ONE_MINUS_COSZERO = 1.0E-10;
	float temp, rnd;


	// DEBUGGING:
	float temp_weight = 0.0f;
	float temp_absorbed = 0.0f;




	const int NUMSTEPS   = 5000;   // Max number of steps photon can take.

	//if (index >= 0 && index <= g_work_size)
	for (i = 0; i < NUMPHOTONS; i++) {         // Photons launched by this work-item (i.e. thread).

		p.x = 0.0f;
		p.y = 0.0f;
		p.z = 0.0f;

		p.status = ALIVE;
		p.weight = 1.0f;
		temp_weight = 1.0f;

		// Randomly set photon trajectory to yield an isotropic source.
		cost = 2.0f*rand(&z1, &z2, &z3, &z4) - 1.0f;
		sint = native_sqrt((float)(1.0f - cost*cost));	// sintheta is always positive
		//sint = sqrt((float)(1.0f - cost*cost));
		psi = 2.0f*PI*rand(&z1, &z2, &z3, &z4);
		p.ux = sint*native_cos((float)psi);
		p.uy = sint*native_sin((float)psi);
		p.uz = cost;

		//while (p.status == ALIVE) {
		for (j = 0; j < NUMSTEPS; j++) {     // Iterations for each photon (i.e. while(ALIVE))
			//============================  HOP  ===================================
			// Step size to take.
			rnd = rand(&z1, &z2, &z3, &z4);
			p.step_size = -native_log(rnd) / (t.mu_a + t.mu_s);

			// Update position of photon.
			p.x += p.step_size * p.ux;
			p.y += p.step_size * p.uy;
			p.z += p.step_size * p.uz;



			// Only perform these steps when scattering has occurred (i.e. step_size > 0)
			//if (p.step_size > 0.0f) {
				//============================  DROP  ===================================
				absorbed = p.weight * (1 - t.albedo);
				temp_absorbed = temp_weight * (1 - t.albedo);
				p.weight -= absorbed;
				temp_weight = temp_weight - temp_absorbed;

				// Using a planar detection geometry.
				r = fabs(p.z);
				ir  = (r/dr);
				if (ir >= NR) ir = NR;
				global_results[ir + (DETECTOR_SIZE * index)] += absorbed;
				//global_results[ir] += absorbed;
				barrier(CLK_LOCAL_MEM_FENCE);



				//============================  SPIN  ===================================
				tempdir = p.ux;

				rnd = rand(&z1, &z2, &z3, &z4);

				if(g==0.0f)
				{
					cost = 2.0f*rnd - 1.0f;              //Should be close close??!!!!!
				}
				else
				{
					float temp1 = (1.0f - g*g) / (1.0f - g + 2*g*rnd);
					cost = (1.0f + g*g - temp1*temp1) / (2.0f*g);
				}
				sint = native_sqrt(1.0f - cost*cost);  // sqrt() is faster than sin().






				// Sample psi.
				psi = 2.0f * PI * rand(&z1, &z2, &z3, &z4);
				cosp = native_cos(psi);
				if (psi < PI)
				{
					sinp = native_sqrt((float)(1.0f - cosp*cosp));
				}
				else
				{
					sinp = -native_sqrt((float)(1.0f - cosp*cosp));
				}




				// New trajectory.
				if ((1 - fabs(p.uz)) <= ONE_MINUS_COSZERO)
				{   // Close to perpendicular.


					uxx = sint * cosp;
					uyy = sint * sinp;
					uzz = cost * (p.uz>=0 ? 1:-1);   //copysign(cost, p.uz*cost);

					temp = native_rsqrt((float)(uxx*uxx + uyy*uyy + uzz*uzz));

					p.ux = uxx * temp;
					p.uy = uyy * temp;
					p.uz = uzz * temp;



				}
				else
				{

					temp = native_sqrt((float)(1.0f - p.uz*p.uz));
					uxx = ((sint * (p.ux * p.uz * cosp - p.uy * sinp)) / temp) + p.ux * cost;
					//p.ux = sint * (p.ux * p.uz * cosp - p.uy * sinp) / temp + p.ux * cost;
					uyy = sint * (p.uy * p.uz * cosp + p.ux * sinp) / temp + p.uy * cost;
					//p.uy = sint * (p.uy * p.uz * cosp + p.ux * sinp) / temp + p.uy * cost;
					uzz = -sint * cosp * temp + p.uz * cost;
					//p.uz = -sint * cosp * temp + p.uz * cost;

					temp = native_rsqrt((float)(uxx*uxx + uyy*uyy + uzz*uzz));

					p.ux = uxx * temp;
					p.uy = uyy*(temp);
					p.uz = uzz*(temp);


				}
				barrier(CLK_LOCAL_MEM_FENCE);




				//============================  ROULETTE  ===================================
				if (p.weight < THRESHOLD)
				{
					if (rand(&z1, &z2, &z3, &z4) <= CHANCE)
					{
						p.weight /= CHANCE;
					}
					else
					{
						p.status = DEAD;
						j = NUMSTEPS;
					}
				}



				barrier(CLK_LOCAL_MEM_FENCE);
			//}



			//============================  ROULETTE  ===================================
			if (p.weight < THRESHOLD)
			{
				if (rand(&z1, &z2, &z3, &z4) <= CHANCE)
				{
					p.weight /= CHANCE;
				}
				else
				{
					p.status = DEAD;
					j = NUMSTEPS;
				}
			}

			barrier(CLK_LOCAL_MEM_FENCE);

		} // end NUMSTEPS loop


		barrier(CLK_LOCAL_MEM_FENCE);


	} // end NUMPHOTONS loop


};





