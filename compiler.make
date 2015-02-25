CXX = g++#-3.4.3#-4.0.2
DEBUG = #-pg #-g #-ggdb
WARNING = -Wall -W -Wno-unused-parameter #-w
OPTIMIZER = -O3#-finline-limit=1000

CXXFLAGS = $(DEBUG) $(WARNING) $(OPTIMIZER)
LIBS = #/usr/lib/libpapi.a
LINK = -lboost_regex -lboost_iostreams
