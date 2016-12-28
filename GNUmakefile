CROSS_COMPILE =
CXX = $(CROSS_COMPILE)g++
CXXFLAGS = -std=gnu++11 -O2 -fPIC -fvisibility=hidden -fno-rtti -fno-exceptions
LDFLAGS =

DESTDIR =
PREFIX = /usr/local

all: libzbessel.so
.PHONY: all

clean:
	rm -f libzbessel.so zbessel.o
.PHONY: clean

install: all
	install -D -m755 libzbessel.so $(DESTDIR)$(PREFIX)/lib/libzbessel.so.1.0
	ln -sf libzbessel.so.1.0 $(DESTDIR)$(PREFIX)/lib/libzbessel.so.1
	ln -sf libzbessel.so.1 $(DESTDIR)$(PREFIX)/lib/libzbessel.so
	install -D -m644 zbessel.h $(DESTDIR)$(PREFIX)/include/zbessel.h
	install -D -m644 zbessel.hh $(DESTDIR)$(PREFIX)/include/zbessel.hh
	cp -dr zbessel $(DESTDIR)$(PREFIX)/include/

libzbessel.so: zbessel.o
	$(CXX) -shared -o $@ $^ $(LDFLAGS) -Wl,-soname,libzbessel.so.1

zbessel.o: zbessel.cc
	$(CXX) $(CXXFLAGS) -c -o $@ $^
