CC       = gcc
CFLAGS   = -std=c11 -Wall -Wextra -Wno-unused-function -O2
CPPFLAGS = -I/usr/include/opencv4
LDLIBS   = -lopencv_core -lopencv_imgproc -lopencv_highgui \
           -lavformat -lavcodec -lswscale -lavutil -lm

TARGET   = imgprc
SRC      = main.c

.PHONY: all clean run

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ $< $(LDLIBS)

run: $(TARGET)
	./$(TARGET) raw/AtomicBlonde.mp4

clean:
	rm -f $(TARGET)
