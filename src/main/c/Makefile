ifeq ($(JAVA_HOME),)
$(error $$(JAVA_HOME) is not defined)
endif

JNI_INCLUDE_DIRS=$(addprefix -I,$(shell find $(JAVA_HOME)/include -type d))
CFLAGS=-O2 -Wall -fPIC $(JNI_INCLUDE_DIRS)
CC=gcc

#OS-dependent extension lookup
UNAME := $(shell uname)
ifeq ($(UNAME),Darwin)
LIB_EXT=Darwin.dylib
else
LIB_EXT=Linux.so
endif

MINIMAP2_COMMIT=v2.17
JNI_BASE_NAME=org_broadinstitute_hellbender_utils_minimap2_MiniMap2Index

all: libmm2.$(LIB_EXT)

libmm2.$(LIB_EXT): $(JNI_BASE_NAME).o minimap2/libminimap2.a
	$(CC) -dynamiclib -shared -o $@ $^ -lm -lz -lpthread

minimap2/libminimap2.a: Makefile
	rm -rf minimap2 && git clone https://github.com/lh3/minimap2 -b $(MINIMAP2_COMMIT) && $(MAKE) CFLAGS="$(CFLAGS)" -C minimap2 libminimap2.a

$(JNI_BASE_NAME).o: $(JNI_BASE_NAME).c version.h minimap2/libminimap2.a

version.h: Makefile
	echo '#define MINIMAP2_VERSION "'$(MINIMAP2_COMMIT)'"' > version.h

clean:
	rm -rf minimap2 *.o version.h libmm2.$(LIB_EXT)

.PHONY: all clean
