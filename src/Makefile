# NIM Makefile

include macros.make

SRCDIR      = $(CURDIR)

all: info build

info:
	$(info Inside top level Makefile)
	$(info ARCH   = $(ARCH))
	$(info OPENMP = $(OPENMP))
	$(info NZ     = $(NZ))
	$(info HW     = $(HW))
	$(info SRCDIR = $(SRCDIR))

build:
	$(MAKE) $(GMAKEMINUSJ) SRCDIR=$(SRCDIR) -C kernel
	$(MAKE) $(GMAKEMINUSJ) SRCDIR=$(SRCDIR) -C diffstats

clean:
	$(MAKE) -C kernel clean
	$(RM) bin/* lib/*


