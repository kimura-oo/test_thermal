
CC = g++
CFLAGS = -Wall -Wextra -pedantic -O3
LDFLAGS =
INCLUDES = -I/usr/local/include -I../../include
LIBS = -L/usr/local/lib -lm -L../../lib -lBB -lBBFE_std -lBBFE_sys -lBBFE_elemmat

TARGET1 = surf_dbc_all
OBJS1 = surf_dbc_all.o surf_core.o
TARGET2 = surf_conn
OBJS2 = surf_conn.o surf_core.o
TARGET3 = surf_nbc
OBJS3 = surf_nbc.o surf_core.o

.SUFFIXES: .c .cpp .o

all: $(TARGET1) $(TARGET2) $(TARGET3)

$(TARGET1): $(OBJS1)
	$(CC) $(CFLAGS) $(INCLUDES) $(LDFLAGS) -o $@ $^ $(LIBS) 

$(TARGET2): $(OBJS2)
	$(CC) $(CFLAGS) $(INCLUDES) $(LDFLAGS) -o $@ $^ $(LIBS) 

$(TARGET3): $(OBJS3)
	$(CC) $(CFLAGS) $(INCLUDES) $(LDFLAGS) -o $@ $^ $(LIBS) 

clean:
		-rm -f $(OBJS1) $(TARGET1) $(OBJS2) $(TARGET2) $(OBJS3) $(TARGET3)

.cpp.o:
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $^ $(LIBS)

.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $^ $(LIBS) 
