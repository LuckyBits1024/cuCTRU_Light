## environment
CXX := g++
CUDAPATH := /usr/local/cuda
CUDA_GENCODE_FLAG := -arch=compute_86 -code=sm_86

NVCC := $(CUDAPATH)/bin/nvcc
CUDAINCLUDEDIR := $(CUDAPATH)/include
CUDALIBDIR := $(CUDAPATH)/lib64

NVCCFLAGS := -std=c++17 --compiler-bindir=$(CXX) $(CUDA_GENCODE_FLAG) -rdc=true \
 -lcuda -lnvToolsExt $(NVCC_DEFINES) -prec-div=false\
-O3 
#-Xnvlink --suppress-stack-size-warning
#--nvlink-options -suppress-stack-size-warning

## sourrce file
## 编译目录放在根目录下
SRCDIR := .
BLDDIR := .

SRCDIR_SOURCE := $(wildcard $(SRCDIR)/*.cu)

FIPS202_WS_DIR := $(SRCDIR)/fips202_ws
FIPS202_WS_SOURCE := $(wildcard $(FIPS202_WS_DIR)/*.cu)

POLY_NTT_653_DIR := $(SRCDIR)/poly_mul_n653q
POLY_NTT_653_SOURCE := $(wildcard $(POLY_NTT_653_DIR)/*.cu)

SOURCE = $(FIPS202_WS_SOURCE) $(SRCDIR_SOURCE) $(POLY_NTT_653_SOURCE)

## command
.PHONY: all clean build

build: $(BLDDIR)/test_ntt

$(BLDDIR)/test_ntt: $(SOURCE)
	$(NVCC) $(NVCCFLAGS) -DFPTRU_N=653 -DNEW=1 -g $^ -o $@ 

$(BLDDIR)/test653_v1: $(SOURCE)
	$(NVCC) $(NVCCFLAGS) -DFPTRU_N=653 -DNEW=0 -g $^ -o $@ 

$(BLDDIR)/test653: $(SOURCE)
	$(NVCC) $(NVCCFLAGS) -DFPTRU_N=653 -DNEW=1 -g $^ -o $@ 

$(BLDDIR)/test761: $(SOURCE)
	$(NVCC) $(NVCCFLAGS) -DFPTRU_N=761 -g $^ -o $@

$(BLDDIR)/test1277: $(SOURCE)
	$(NVCC) $(NVCCFLAGS) -DFPTRU_N=1277 -g $^ -o $@


list = 1 2 4 8 16 32 64 128 256 512 1024
batch_thread: $(SOURCE)
	for i in $(list); do \
		for j in $$(seq 1 10); do \
			$(NVCC) $(NVCCFLAGS) -DFPTRU_N=653 -DNEW=1 -DBATCH_SIZE=$$i -DNUM_THREAD=$$j -g $^ -o ./hxw/batch_$${i}_thread_$${j}; \
		done; \
	done


clean:
#	rm -r $(BLDDIR)/test1277
#	rm -r $(BLDDIR)/test761
#	rm -r $(BLDDIR)/test653_v1
	rm -r $(BLDDIR)/test_ntt
#	rm -r $(BLDDIR)/test653


