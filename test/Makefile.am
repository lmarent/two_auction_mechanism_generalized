check_PROGRAMS = test_runner

API_INC			= $(top_srcdir)/include
INC_DIR 		= $(top_srcdir)/include/


test_runner_SOURCES =  @top_srcdir@/src/Error.cpp \
					   @top_srcdir@/src/TwoAuctionMechanismGeneralized.cpp \
					   @top_srcdir@/test/TwoAuctionMechanismGeneralized_test.cpp \
					   @top_srcdir@/test/test_runner.cpp

test_runner_CPPFLAGS  = -I$(API_INC) $(CPPUNIT_CFLAGS) 

test_runner_LDADD  = $(CPPUNIT_LIBS) -ldl -lcppunit 

TESTS = $(check_PROGRAMS)

if ENABLE_DEBUG
  AM_CXXFLAGS = -I$(top_srcdir)/include \
				-g  -fno-inline -DDEBUG 
else
  AM_CXXFLAGS = -I$(top_srcdir)/include \
				-O2 
endif

AM_CXXFLAGS += -Wall -ggdb -pedantic -Wno-long-long


#  TEST_CXX_FLAGS=`cppunit-config --cflags`
