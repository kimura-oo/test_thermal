
CC = mpic++
CFLAGS = -Wall -Wextra -pedantic -O3
LDFLAGS =
INCLUDES = -I/usr/local/include -I../include -I../submodule/monolis/include
LIBS = -L/usr/local/lib -lm -L../lib/BB -lBB -lBBFE_std -lBBFE_sys -lBBFE_elemmat

TARGET_LIB = libBBFE_manusol.a
#OBJS = debug_main.o
OBJS = manusol.o

$(TARGET_LIB): $(OBJS)
		ar rcs $@ $?

clean:
		-rm -f $(OBJS) $(TARGET_LIB)

.cpp.o:
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $^ $(INCLUDES)

.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $^ $(INCLUDES)
