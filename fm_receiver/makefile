CFLAGS=-std=c99 -g -W `pkg-config --cflags librtlsdr`
LDLIBS=`pkg-config --libs librtlsdr` -lfftw3 -lm
all: fm_receiver

fm_receiver:

test: fm_receiver
	./fm_receiver -f 91500000 > foo.raw
	sox -r 44100 -e signed -b 16 -c 1 -L foo.raw foo.wav
	mplayer foo.wav
