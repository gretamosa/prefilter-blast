# MAKEFILE NO VALE

CC = gcc
CFLAGS = -Wall -g
LFLAGS = -lm
RUTAUTIL = util

NVCC = /opt/cuda/bin/nvcc
LINK = g++
CUDAPATH = /usr/local/cuda
SDKPATH = /home/telefonica/NVIDIA_GPU_Computing_SDK/C
CUDALIBS = -L$(CUDAPATH)/lib64 -L$(SDKPATH)/lib -L$(SDKPATH)/common/lib/linux -L/home/telefonica/NVIDIA_GPU_Computing_SDK/shared/lib -lcuda -lcudart -lcutil_x86_64 -lshrutil_x86_64

all: Serie

Serie: prefilterSerie.o parse.o
	$(CC) $(CFLAGS) $(LFLAGS) prefilterSerie.o parse.o -o prefS

Paralel: prefilterCuda.o prefilterParalel.o parse.o
	$(LINK) prefilterCuda.o prefilterParalel.o parse.o -o prefP $(CUDALIBS)

util: parse.o

# util
parse.o: $(RUTAUTIL)/parse.c $(RUTAUTIL)/parse.h
	$(CC) $(CFLAGS) $(LFLAGS) -c $(RUTAUTIL)/parse.c

#prefilters
prefilterParalel.o: prefilterParalel.c
	$(CC) $(CFLAGS) $(LFLAGS) -c prefilterParalel.c
	
prefilterSerie.o: prefilterSerie.c
	$(CC) $(CFLAGS) $(LFLAGS) -c prefilterSerie.c	
	
#cuda
prefilterCuda.o: prefilterCuda.cu 
	$(NVCC) -c -g prefilterCuda.cu

clean:
	rm *.o prefS prefP
