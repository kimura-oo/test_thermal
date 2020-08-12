CC = mpic++
FC = mpif90

CFLAGS = -g -O3
LDFLAGS =
INCLUDES = -I/usr/local/include -I./submodule/monolis/include
LIBS = -L/usr/local/lib -lm -lstdc++ -L./submodule/monolis/lib -lmonolis -lmetis

TARGET1 = convdiff_s_msol
TARGET2 = convdiff_ns_msol
OBJS1 = convdiff_s_msol.o convdiff_core.o
OBJS2 = convdiff_ns_msol.o convdiff_core.o

LIBOBJS = libBB/std.o libBB/calc.o libBB/vtk.o \
	   FE_std/integ.o FE_std/shapefunc.o FE_std/mapping.o FE_std/surface.o \
	   FE_sys/memory.o FE_sys/read.o FE_sys/write.o FE_sys/monowrap.o\
	   FE_elemmat/set.o FE_elemmat/equivval.o FE_elemmat/convdiff.o\
	   FE_manusol/manusol.o

.SUFFIXES: .c .cpp .o

all: $(TARGET1) $(TARGET2)

$(TARGET1): $(OBJS1) $(LIBOBJS)
	$(FC) $(CFLAGS) $(INCLUDES) $(LDFLAGS) -o $@ $^ $(LIBS)

$(TARGET2): $(OBJS2) $(LIBOBJS)
	$(FC) $(CFLAGS) $(INCLUDES) $(LDFLAGS) -o $@ $^ $(LIBS)

clean:
		-rm -f $(LIBOBJS) $(OBJS1) $(TARGET1) $(OBJS2) $(TARGET2)

.cpp.o:
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $^

.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $^
