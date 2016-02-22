#
# FLAGS
#

PROGRAM_NAME = cta
VERSION_MAJOR = 0
VERSION_MINOR = 1
VERSION_PATCH = 0

VERSION = $(VERSION_MAJOR).$(VERSION_MINOR).$(VERSION_PATCH)

# PATHS
SRC_DIR = src
BUILD_DIR = build

SRC_CPP := $(wildcard $(SRC_DIR)/*.cpp)
OBJECTS := $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(SRC_CPP))

CXXFLAGS = -std=c++11 -I.
LDLIBS = -lboost_iostreams -lz
LDFLAGS = -static
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
	CXXFLAGS += -D LINUX -pthread  $(INCLUDES) -O3 -g
	ifdef BOOST_MODULE
		INCLUDES += -I $(BOOST_MODULE)/include
		LDFLAGS += -L $(BOOST_MODULE)/lib
	endif
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

PREFIX = /usr/local

#
# TARGETS
#

all: checkdirs $(BUILD_DIR)/cta

install: $(BUILD_DIR)/cta
	install -m 0755 $(BUILD_DIR)/cta $(PREFIX)/bin

clean:
	rm -rf $(BUILD_DIR)

checkdirs: $(BUILD_DIR)

$(BUILD_DIR):
	@mkdir -p $@

$(BUILD_DIR)/cta: $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $^ $(LDFLAGS) $(LDLIBS)

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
