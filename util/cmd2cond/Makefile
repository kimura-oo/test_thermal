
CC = g++
CFLAGS = -Wall -Wextra -pedantic -O3
LDFLAGS =
INCLUDES = -I/usr/local/include -I../../include
LIBS = -L/usr/local/lib -lm -L../../lib -lBB

TARGET1 = cmd2cond
OBJS1 = cmd2cond.o

.SUFFIXES: .c .cpp .o

all: $(TARGET1) $(TARGET2)

$(TARGET1): $(OBJS1)
	$(CC) $(CFLAGS) $(INCLUDES) $(LDFLAGS) -o $@ $^ $(LIBS) 

clean:
		-rm -f $(OBJS1) $(TARGET1) $(OBJS2) $(TARGET2)

.cpp.o:
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $^ $(LIBS)

.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $^ $(LIBS) 
