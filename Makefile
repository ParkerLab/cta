#
# FLAGS
#


CXXFLAGS = -std=c++11 -I.
LIBS = -lboost_iostreams -lz
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)  # Lab server: RHEL6
	INCLUDES += -I $(BOOST_MODULE)/include
	CXXFLAGS += -D LINUX -pthread  $(INCLUDES) -O3 -g
	LDFLAGS = -L $(BOOST_MODULE)/lib $(LIBS) # -static
endif
ifeq ($(UNAME_S),Darwin)  # Mac with Homebrew
	CXXFLAGS += -D OSX -O3
	LDFLAGS = $(LIBS)
endif

UNAME_P := $(shell uname -p)
ifeq ($(UNAME_P),x86_64)
	CXXFLAGS += -D AMD64
endif
ifneq ($(filter %86,$(UNAME_P)),)
	CXXFLAGS += -D IA32
endif
ifneq ($(filter arm%,$(UNAME_P)),)
	CXXFLAGS += -D ARM
endif

PREFIX = /usr/local

#
# TARGETS
#

all: cta

clean:
	rm -f cta

cta: cta.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $^ $(LDFLAGS)

install:
	install -m 0755 cta $(PREFIX)/bin
