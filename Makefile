#
# FLAGS
#

PROGRAM_NAME = cta
VERSION_MAJOR = 0
VERSION_MINOR = 1
VERSION_PATCH = 3

VERSION = $(VERSION_MAJOR).$(VERSION_MINOR).$(VERSION_PATCH)

# PATHS
SRC_DIR = src
BUILD_DIR = build

SRC_CPP := $(wildcard $(SRC_DIR)/*.cpp)
OBJECTS := $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(SRC_CPP))

CPPFLAGS = -pedantic -Wall -Wextra -Wwrite-strings -Wstrict-overflow -Wp,-D_FORTIFY_SOURCE=2 -fexceptions -fstack-protector-strong -fno-strict-aliasing -fPIC $(INCLUDES)
CXXFLAGS = -std=c++11 -O3 -g $(CPPFLAGS)

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
	CXXFLAGS += -D LINUX $(INCLUDES) -O3 -g
endif
ifeq ($(UNAME_S),Darwin)  # Mac with Homebrew
	CXXFLAGS += -D OSX -O3
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

ifdef BOOST_TAGGED
	BOOST_LIBS = -lboost_filesystem-mt -lboost_iostreams-mt -lboost_system-mt -lboost_chrono-mt -lz
else
	BOOST_LIBS = -lboost_filesystem -lboost_iostreams -lboost_system -lboost_chrono -lz
endif

ifdef BOOST_ROOT
	CPPFLAGS += -I$(BOOST_ROOT)/include
	LDFLAGS += -L$(BOOST_ROOT)/lib
else
	ifdef BOOST_INCLUDE
		CPPFLAGS += -I$(BOOST_INCLUDE
	endif

	ifdef BOOST_LIB
		LDFLAGS += -L$(BOOST_LIB)
	endif
endif

LDLIBS = $(BOOST_LIBS)

PREFIX = /usr/local

#
# TARGETS
#

.PHONY: all checkdirs clean install

all: $(BUILD_DIR)/cta

checkdirs: $(BUILD_DIR)

clean:
	rm -rf $(BUILD_DIR)

install: $(BUILD_DIR)/cta
	install -d -m 0755 $(PREFIX)/bin
	install -m 0755 $(BUILD_DIR)/cta $(PREFIX)/bin/cta

$(BUILD_DIR):
	@mkdir -p $@

$(BUILD_DIR)/cta: checkdirs $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $(OBJECTS) $(LDFLAGS) $(LDLIBS)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp $(SRC_DIR)/Version.hpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ -c $<

$(SRC_DIR)/Version.hpp: Makefile
	@echo '//' >$@
	@echo '// Copyright 2015 The Parker Lab at the University of Michigan' >>$@
	@echo '//' >>$@
	@echo '// Licensed under Version 3 of the GPL or any later version' >>$@
	@echo '//' >>$@
	@echo >>$@
	@echo '#define PROGRAM_NAME "$(PROGRAM_NAME)"' >>$@
	@echo '#define VERSION_MAJOR $(VERSION_MAJOR)' >>$@
	@echo '#define VERSION_MINOR $(VERSION_MINOR)' >>$@
	@echo '#define VERSION_PATCH $(VERSION_PATCH)' >>$@
