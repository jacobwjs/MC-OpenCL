

//
//

////////////////////////////////////////////////////////////////////////////////

#include <assert.h>
#include <math.h>
#include <time.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <CL/cl.h>
////////////////////////////////////////////////////////////////////////////////

// Use a static data size for simplicity
//

////////////////////////////////////////////////////////////////////////////////


#define DATA_SIZE (101)
#define MAX_THREADS (128)
#define LOCAL_SIZE (128)			// The local size of a work group.
#define NUM_EVENTS (1)			// Number of times the kernel is executed.
#define PHOTONS_PER_ITEM (4000)	// Number of photons to simulate per-work-item.

double radial_size = 3.0;
int NR = 100;

char * load_program_source(const char *filename);
void printFluence(float *results, const int numPhotons);
void printAllResults(float *results);


int main(int argc, char** argv)
{
	cl_int err;                         // error code returned from api calls

	float results[DATA_SIZE * MAX_THREADS];           // results returned from device

	size_t global;                      // global domain size for our calculation
	size_t local;                       // local domain size for our calculation


	cl_context 		context;            // compute context
	cl_command_queue commands;          // compute command queue
	cl_program 		program;            // compute program
	cl_kernel 		kernel;             // compute kernel

	cl_mem 			input;              // device memory used for the input array
	cl_mem 			output;             // device memory used for the output array
	cl_mem			scatt_coeffs;		// holds various scattering coefficients for simulation.

	cl_event event[NUM_EVENTS];			// Events executed on the device.

	float tissue_scatt_coeffs[MAX_THREADS];


	//Get an OpenCL platform
	cl_platform_id 	cpPlatform;
	clGetPlatformIDs(1, &cpPlatform, NULL);

	// Get a GPU device
	cl_device_id cdDevice;
	err = clGetDeviceIDs(cpPlatform, CL_DEVICE_TYPE_GPU, 1, &cdDevice, NULL);
	if (err != CL_SUCCESS)
	{
		printf("Error: Failed to create a device group!\n");
		return EXIT_FAILURE;
	}


	// Holds the random numbers used for seeding the RNG in the kernel.  Each work-item
	// produces it's own independent stream of random numbers, and since the RNG being
	// used needs 4 seeds, we generate 4 * workitems random seeds using the RNG library
	// from the C-library.
	//srand(time(0));
	srand(11);
	unsigned int rng_seeds[4*MAX_THREADS];
	int i;
	for(i = 0; i < 4*MAX_THREADS; i++)
		rng_seeds[i] = rand() + 128 + rand();


	// Initialize results to zero.
	for (i = 0; i < DATA_SIZE; i++)
		results[i] = 0;

	// Initialize the scattering coefficients of the tissue.
	for (i = 0; i < MAX_THREADS; i++) {
		//tissue_scatt_coeffs[i] = rand() % 100;
		tissue_scatt_coeffs[i] = 100.0f + 0;//(i*0.01f);
		printf("scatter coeffs = %f\n", tissue_scatt_coeffs[i]);
	}
	tissue_scatt_coeffs[0] = 100.0f;


	char cBuffer[1024];
	clGetDeviceInfo(cdDevice, CL_DEVICE_NAME, sizeof(cBuffer), &cBuffer, NULL);
	printf("CL_DEVICE_NAME:       %s\n", cBuffer);
	clGetDeviceInfo(cdDevice, CL_DRIVER_VERSION, sizeof(cBuffer), &cBuffer, NULL);
	printf("CL_DRIVER_VERSION: %s\n\n", cBuffer);



	// Create a compute context
	//
	context = clCreateContext(0, 1, &cdDevice, NULL, NULL, &err);
	if (!context)
	{
		printf("Error: Failed to create a compute context!\n");
		return EXIT_FAILURE;
	}

	// Create a command commands
	//
	commands = clCreateCommandQueue(context, cdDevice, CL_QUEUE_PROFILING_ENABLE, &err);
	if (!commands)
	{
		printf("Error: Failed to create a command commands!\n");
		return EXIT_FAILURE;
	}

	// Create the compute program from the source buffer
	//
	const char *filename = "launchPhotons.opencl";
	char *program_source = load_program_source(filename);
	program = clCreateProgramWithSource(context, 1, (const char **)&program_source, NULL, &err);
	if (!program)
	{
		printf("Error: Failed to create compute program!\n");
		return EXIT_FAILURE;
	}

	// Build the program executable
	//
	err = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);
	if (err != CL_SUCCESS)
	{
		size_t len;
		char buffer[2048];

		printf("Error: Failed to build program executable!\n");
		clGetProgramBuildInfo(program, cdDevice, CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);
		printf("%s\n", buffer);
		exit(1);
	}

	// Create the compute kernel in the program we wish to run
	//
	kernel = clCreateKernel(program, "LaunchPhoton", NULL);
	if (!kernel || err != CL_SUCCESS)
	{
		printf("Error: Failed to create compute kernel!\n");
		exit(1);
	}

	// Create the input and output arrays in device memory for our calculation
	//
	input  = clCreateBuffer(context,  CL_MEM_READ_ONLY,  sizeof(unsigned int) * (MAX_THREADS*4), NULL, NULL);
	output = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(float) * DATA_SIZE * MAX_THREADS, NULL, NULL);
	scatt_coeffs  = clCreateBuffer(context,  CL_MEM_READ_ONLY,  sizeof(float) * MAX_THREADS, NULL, NULL);
	if (!input || !output)
	{
		printf("Error: Failed to allocate device memory!\n");
		exit(1);
	}

	// Write our data set into the input array in device memory
	//
	err = clEnqueueWriteBuffer(commands, input, CL_TRUE, 0, sizeof(unsigned int) * (MAX_THREADS*4), rng_seeds, 0, NULL, NULL);
	err |= clEnqueueWriteBuffer(commands, scatt_coeffs, CL_TRUE, 0, sizeof(float) * MAX_THREADS, tissue_scatt_coeffs, 0, NULL, NULL);
	if (err != CL_SUCCESS)
	{
		printf("Error: Failed to write to source array!\n");
		exit(1);
	}

	// Set the arguments to our compute kernel
	//
	err = 0;
	const int num_photons_per_work_item = PHOTONS_PER_ITEM; ///MAX_THREADS;
	const int detector_size = DATA_SIZE;
	//err  = clSetKernelArg(kernel, 0, sizeof(cl_mem), &input);
	//err |= clSetKernelArg(kernel, 1, sizeof(cl_mem), &output);
	//err |= clSetKernelArg(kernel, 2, sizeof(unsigned int), &count);
	err  = clSetKernelArg(kernel, 0, sizeof(cl_mem), &input);
	err |= clSetKernelArg(kernel, 1, sizeof(cl_mem), &output);
	err |= clSetKernelArg(kernel, 2, sizeof(const int), &num_photons_per_work_item);
	err |= clSetKernelArg(kernel, 3, sizeof(const int), &detector_size);
	err |= clSetKernelArg(kernel, 4, sizeof(cl_mem), &scatt_coeffs);
	if (err != CL_SUCCESS)
	{
		printf("Error: Failed to set kernel arguments! %d\n", err);
		exit(1);
	}

	// Get the maximum work group size for executing the kernel on the device
	//
	err = clGetKernelWorkGroupInfo(kernel, cdDevice, CL_KERNEL_WORK_GROUP_SIZE, sizeof(local), &local, NULL);
	if (err != CL_SUCCESS)
	{
		printf("Error: Failed to retrieve kernel work group info! %d\n", err);
		exit(1);
	}

	printf("local work-group-size = %d\n", (int)local);

	// Execute the kernel over the entire range of our 1d input data set
	// using the maximum number of work group items for this device
	//
	global = MAX_THREADS;
	local = LOCAL_SIZE;

	for (i = 0; i < NUM_EVENTS; i++) {
		err  = clEnqueueNDRangeKernel(commands, kernel, 1, NULL, &global, &local, 0, NULL, &event[i]);
	}

	if (err)
	{
		printf("Error: Failed to execute kernel!\n");
		return EXIT_FAILURE;
	}
	// Wait for the command commands to get serviced before reading back results
	//
	//clFinish(commands);
	clWaitForEvents((cl_uint)i, event);


	// Read back the results from the device to verify the output
	//
	err = clEnqueueReadBuffer( commands, output, CL_TRUE, 0, sizeof(float) * DATA_SIZE * MAX_THREADS, results, 0, NULL, NULL );
	if (err != CL_SUCCESS)
	{
		printf("Error: Failed to read output array! %d\n", err);
		exit(1);
	}

	// Wait for the command commands to get serviced before reading back results
	//
	clFinish(commands);

	cl_ulong start, end;
	clGetEventProfilingInfo(event[0], CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, NULL);
	clGetEventProfilingInfo(event[i-1], CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &end, NULL);
	printf("Time for event (ms): %10.5f  \n", (end-start)/1000000.0);


	// Display the total calculations made on the GPU.
	printAllResults(results);

	// Sum up all values.
	float vals = 0.0f;
	float summed_results[DATA_SIZE];
	int j;
	for (i = 0; i < DATA_SIZE; i++) {
		for (j = 0; j < MAX_THREADS; j++) {
			vals += results[(j*DATA_SIZE) + i];
		}
		summed_results[i] = vals;
		vals = 0.0f;
	}

	// Print out the calculated fluences to disk.
	printFluence((float *)summed_results, (MAX_THREADS)*num_photons_per_work_item);


	// Shutdown and cleanup
	//
	clReleaseMemObject(input);
	clReleaseMemObject(output);
	clReleaseProgram(program);
	clReleaseKernel(kernel);
	clReleaseCommandQueue(commands);
	clReleaseContext(context);

	for (i = 0; i < NUM_EVENTS; i++)
	{
		clReleaseEvent(event[i]);
	}


	return 0;
}



char * load_program_source(const char *filename)
		{
	struct stat statbuf;
	FILE *fp;
	char *source;
	fp = fopen(filename, "r");
	printf("OpenCL file = %s\n", filename);
	assert (fp != NULL);
	stat(filename, &statbuf);
	source = (char *)malloc(statbuf.st_size + 1);
	fread(source, statbuf.st_size, 1, fp);
	source[statbuf.st_size] = '\0';

	//printf("kernel = %s", source);

	return source;
}

void printAllResults(float *results)
{
	FILE *target; 	// File pointer to save data to.
	target = fopen("debug.txt", "w");

	int i, j;
	for (i = 0; i < DATA_SIZE; i++) {
		for (j = 0; j < MAX_THREADS; j++) {
			fprintf(target, "[%d, %d] = %5.3f \t ", j, i, results[(j*DATA_SIZE) + i]);
		}
		fprintf(target,"\n");
	}
	fclose(target);
}



void printFluence(float *results, const int numPhotons)
{
	FILE *target; 	// File pointer to save data to.
	target = fopen("fluence.txt", "w");

	printf("Simulated %d photons...\n", numPhotons);

	int ir;
	double shellvolume;
	double dr = radial_size / NR;	// cm
	double r;
	double Fpla;
	double Nphotons = (numPhotons);///(32);
	double mu_a = 1.0;

	for (ir = 0; ir <= NR; ir++)
	{
		r = (ir + 0.5)*dr;
		shellvolume = dr;	// per cm^2 area of plane.
		Fpla = results[ir]/Nphotons/shellvolume/mu_a;
		printf("%5.8f \t %5.5f \t %4.3e \n", results[ir], r, Fpla);
		fprintf(target, "%5.5f \t %4.3e \n", r, Fpla);
	}

	fclose(target);

}


