CC = g++
CFLAGS = -O3
LDFLAGS =
INCLUDES = -I/usr/local/include
LIBS = -L/usr/local/lib -lm 

TARGET = thermal
OBJS = main.o solve_mat.o

.SUFFIXES: .c .cpp .o

all: $(TARGET)

$(TARGET): $(OBJS)
		$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^ $(LIBS)

clean:
		-rm -f $(OBJS)

.cpp.o:
	$(CC) $(CFLAGS) -o $@ -c $^ 
