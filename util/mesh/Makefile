
CC = g++
CFLAGS = -Wall -Wextra -pedantic -O3
LDFLAGS =
INCLUDES = -I/usr/local/include -I../../include
LIBS = -L/usr/local/lib -lm -L../../lib -lBB -lBBFE_std -lBBFE_sys

TARGET1 = mesh_extract
OBJS1 = mesh_extract.o mesh_core.o

TARGET2 = mesh_remove
OBJS2 = mesh_remove.o mesh_core.o

TARGET3 = mesh_surf_extract
OBJS3 = mesh_surf_extract.o mesh_core.o

TARGET4 = mesh_surf_remove
OBJS4 = mesh_surf_remove.o mesh_core.o

.SUFFIXES: .c .cpp .o

all: $(TARGET1) $(TARGET2) $(TARGET3) $(TARGET4)

$(TARGET1): $(OBJS1)
	$(CC) $(CFLAGS) $(INCLUDES) $(LDFLAGS) -o $@ $^ $(LIBS) 

$(TARGET2): $(OBJS2)
	$(CC) $(CFLAGS) $(INCLUDES) $(LDFLAGS) -o $@ $^ $(LIBS) 

$(TARGET3): $(OBJS3)
	$(CC) $(CFLAGS) $(INCLUDES) $(LDFLAGS) -o $@ $^ $(LIBS) 

$(TARGET4): $(OBJS4)
	$(CC) $(CFLAGS) $(INCLUDES) $(LDFLAGS) -o $@ $^ $(LIBS) 

clean:
		-rm -f $(OBJS1) $(TARGET1) $(OBJS2) $(TARGET2) $(OBJS3) $(TARGET3) $(OBJS4) $(TARGET4)

.cpp.o:
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $^ $(LIBS)

.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $^ $(LIBS) 
