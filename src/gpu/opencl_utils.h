
/**
 * Copyright (c) 2011 MOSAIC Group (ETH Zurich)
 * 
 * Author: Omar Awile
 *
 * Note:
 * This software contains source code provided by NVIDIA Corporation.
 **/
#include <CL/opencl.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>


void tic(struct timeval* t);

double toc(struct timeval* t);

void clReadCode(char* fileName, char** code, size_t* len);

int clInit(cl_platform_id *platform, cl_device_id *device, cl_context *context, cl_command_queue *queue, cl_program *program);

int clBuild(char **filename, int fcount, cl_program *program, cl_context *context, cl_device_id *device, const char* option);

const char* clErrorString(cl_int error);

void checkError(cl_int error, const char* string, const char* funct_name, int line, char* stamp);

int divide_roundUp(int toBeRounded, int divider);
