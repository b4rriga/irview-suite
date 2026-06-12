CC = gcc
CFLAGS = -O3 -march=native -ffast-math -fPIC
LDFLAGS = -shared

TARGET = ircore.so
SRC = ircore.c

PREFIX ?= /usr/local/bin

all: $(TARGET)

$(TARGET):
	$(CC) $(CFLAGS) $(LDFLAGS) -o $(TARGET) $(SRC) -lm

install: all
	install -Dm755 main.py $(PREFIX)/irview
	install -Dm755 ircap.py $(PREFIX)/ircap
	install -Dm755 irshot.py $(PREFIX)/irshot
	install -Dm755 irwebcam.py $(PREFIX)/irwebcam
	install -Dm644 config.ini /etc/irview.ini

uninstall:
	rm -f $(PREFIX)/irview
	rm -f $(PREFIX)/ircap
	rm -f $(PREFIX)/irshot
	rm -f $(PREFIX)/irwebcam
	rm -f /etc/irview.ini

clean:
	rm -f $(TARGET)
