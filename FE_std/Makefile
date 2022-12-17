
CC = mpic++
CFLAGS = -Wall -Wextra -pedantic -O3
LDFLAGS =
INCLUDES = -I/usr/local/include -I../include
LIBS = -L/usr/local/lib -lm -L../lib/BB -lBB

TARGET_LIB = libBBFE_std.a
#OBJS = debug_main.o
OBJS = integ.o shapefunc.o shapefunc_gen.o mapping.o surface.o

$(TARGET_LIB): $(OBJS)
		ar rcs $@ $?

clean:
		-rm -f $(OBJS) $(TARGET_LIB)

.cpp.o:
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $^ $(LIBS)

.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $^ $(LIBS)
