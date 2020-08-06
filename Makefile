CC = mpic++
FC = mpif90

CFLAGS = -g -O3
LDFLAGS =
INCLUDES = -I/usr/local/include -I./submodule/monolis/include
LIBS = -L/usr/local/lib -lm -lstdc++ -L./submodule/monolis/lib -lmonolis -lmetis

TARGET = thermal
OBJS = main.o

.SUFFIXES: .c .cpp .o

all: $(TARGET)

$(TARGET): $(OBJS)
	$(FC) $(CFLAGS) $(INCLUDES) $(LDFLAGS) -o $@ $^ $(LIBS)

clean:
		-rm -f $(OBJS)

.cpp.o:
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $^

.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $^
