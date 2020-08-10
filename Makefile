CC = mpic++
FC = mpif90

CFLAGS = -g -O3
LDFLAGS =
INCLUDES = -I/usr/local/include -I./submodule/monolis/include
LIBS = -L/usr/local/lib -lm -lstdc++ -L./submodule/monolis/lib -lmonolis -lmetis

TARGET = convdiff_s
OBJS = convdiff_s.o \
	   libBB/std.o libBB/calc.o libBB/vtk.o \
	   FE_std/integ.o FE_std/shapefunc.o FE_std/mapping.o FE_std/surface.o \
	   FE_sys/memory.o FE_sys/read.o FE_sys/write.o FE_sys/monowrap.o\
	   FE_elemmat/set.o FE_elemmat/equivval.o FE_elemmat/convdiff.o\
	   FE_manusol/manusol.o

.SUFFIXES: .c .cpp .o

all: $(TARGET)

$(TARGET): $(OBJS)
	$(FC) $(CFLAGS) $(INCLUDES) $(LDFLAGS) -o $@ $^ $(LIBS)

clean:
		-rm -f $(OBJS) $(TARGET)

.cpp.o:
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $^

.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $^
