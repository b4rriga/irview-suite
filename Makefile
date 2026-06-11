CC = gcc
CFLAGS = -O3 -march=native -ffast-math -fPIC
LDFLAGS = -shared

TARGET = ircore.so
SRC = ircore.c

PREFIX ?= /usr/local
BIN_NAME ?= irview
BIN_PATH = $(PREFIX)/bin/$(BIN_NAME)

all: $(TARGET)

$(TARGET):
	$(CC) $(CFLAGS) $(LDFLAGS) -o $(TARGET) $(SRC) -lm

install: all
	install -Dm755 main.py $(BIN_PATH)

uninstall:
	rm -f $(BIN_PATH)

clean:
	rm -f $(TARGET)
