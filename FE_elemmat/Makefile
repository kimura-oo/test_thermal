
CC = mpic++
CFLAGS = -Wall -Wextra -pedantic -O3
LDFLAGS =
INCLUDES = -I/usr/local/include -I../submodule/monolis/include -I../include
LIBS = -L/usr/local/lib -lm -L../lib/BB -lBB -lBBFE_std

TARGET_LIB = libBBFE_elemmat.a
#OBJS = debug_main.o
OBJS = set.o equivval.o convdiff.o fluid.o solid.o

$(TARGET_LIB): $(OBJS)
		ar rcs $@ $?

clean:
		-rm -f $(OBJS) $(TARGET_LIB)

.cpp.o:
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $^ $(LIBS)

.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $^ $(LIBS)
