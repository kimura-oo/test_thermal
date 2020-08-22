CC = g++
CFLAGS = -O3
LDFLAGS =
INCLUDES = -I/usr/local/include
LIBS = -L/usr/local/lib -lm 

TARGET1 = meshgen_tet
OBJS1 = meshgen_tet.o

TARGET2 = meshgen_hex
OBJS2 = meshgen_hex.o

.SUFFIXES: .c .cpp .o

all: $(TARGET1) $(TARGET2)

$(TARGET1): $(OBJS1)
		$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^ $(LIBS)
$(TARGET2): $(OBJS2)
		$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^ $(LIBS)

clean:
		-rm -f $(OBJS1) $(OBJS2)

.cpp.o:
	$(CC) $(CFLAGS) -o $@ -c $^ 
