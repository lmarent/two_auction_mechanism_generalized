# -----------------------------------*- mode: Makefile; -*--
# Makefile.am - Makefile.am for  Two Auction Generalized Mechanism library
# ----------------------------------------------------------
# $Id$
# $HeadURL$
# ==========================================================
#                      
# (C)opyright, all rights reserved by
# - System and Computing Engineering, Universidad de los Andes
# ==========================================================
#
#


pkgincludedir = $(includedir)/$(PACKAGE)

API_INC		= $(top_srcdir)/include
INC_DIR 	= $(top_srcdir)/include

lib_LTLIBRARIES = libtwoauctiongen.la

if ENABLE_DEBUG
AM_CXXFLAGS = -I$(top_srcdir)/include -g  -fno-inline -DDEBUG -ggdb
else
AM_CXXFLAGS = -I$(top_srcdir)/include -O2
endif


#-----------------------------------------------------------------------
#--------------------------- Two Auction Generalized 
#-----------------------------------------------------------------------

libtwoauctiongen_la_CPPFLAGS = 

libtwoauctiongen_la_SOURCES = Error.cp \
							  TwoAuctionMechanismGeneralized.cpp \
							  TwoAuctionMechanismGeneralizedWrapper.cpp 


# end of Makefile.am

