
CC = g++
CFLAGS = -Wall -Wextra -pedantic -O3
LDFLAGS =
INCLUDES = -I/usr/local/include -I../../include 
LIBS = -L/usr/local/lib -lm -L../../lib -lBB -lBBFE_std -lBBFE_sys -lBBFE_elemmat

TARGET1 = surf_dbc_all
OBJS1 = surf_dbc_all.o surf_core.o
TARGET2 = surf_conn
OBJS2 = surf_conn.o surf_core.o
TARGET3 = surf_dbc
OBJS3 = surf_dbc.o surf_core.o
TARGET4 = surf_nbc
OBJS4 = surf_nbc.o surf_core.o
TARGET5 = surf_bc_merge
OBJS5 = surf_bc_merge.o surf_core.o

.SUFFIXES: .c .cpp .o

all: $(TARGET1) $(TARGET2) $(TARGET3) $(TARGET4) $(TARGET5)

$(TARGET1): $(OBJS1)
	$(CC) $(CFLAGS) $(INCLUDES) $(LDFLAGS) -o $@ $^ $(LIBS) 

$(TARGET2): $(OBJS2)
	$(CC) $(CFLAGS) $(INCLUDES) $(LDFLAGS) -o $@ $^ $(LIBS) 

$(TARGET3): $(OBJS3)
	$(CC) $(CFLAGS) $(INCLUDES) $(LDFLAGS) -o $@ $^ $(LIBS) 

$(TARGET4): $(OBJS4)
	$(CC) $(CFLAGS) $(INCLUDES) $(LDFLAGS) -o $@ $^ $(LIBS) 

$(TARGET5): $(OBJS5)
	$(CC) $(CFLAGS) $(INCLUDES) $(LDFLAGS) -o $@ $^ $(LIBS) 

clean:
		-rm -f $(OBJS1) $(TARGET1) $(OBJS2) $(TARGET2) $(OBJS3) $(TARGET3) $(OBJS4) $(TARGET4) $(OBJS5) $(TARGET5)

.cpp.o:
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $^ $(LIBS)

.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $^ $(LIBS) 
