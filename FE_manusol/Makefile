
CC = mpic++
CFLAGS = -Wall -Wextra -pedantic -O3
LDFLAGS =
INCLUDES = -I/usr/local/include
LIBS = -L/usr/local/lib -lm

TARGET_LIB = ../lib/libBBFE_manusol.a
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
