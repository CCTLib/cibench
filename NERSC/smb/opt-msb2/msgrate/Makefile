C = mpicc -cc=gcc
#
# Copyright 2006 Sandia Corporation. Under the terms of Contract
# DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
# retains certain rights in this software.
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, 
# Boston, MA  02110-1301, USA.

RELEASE_VERSION = 1.0

# Please have a valid MPI compiler installed and set following path
MPICC = 
CC = mpicc -cc=gcc -g -O3 -flto
#CC1 = mpicc -cc=gcc -g -O3 -flto -fprofile-generate
#CC2 = mpicc -cc=gcc -g -O3 -flto -fprofile-use -fprofile-correction
MPI_COMPILE_FLAGS = -I$(MPICC_INCLUDE)
MPI_LINK_FLAGS = -L$(MPICC_LIBRARY) -Wl,-rpath -Wl,$(MPICC_LABRARY) -Wl,--enable-new-dtags -lmpi

ifndef MPICC
$(error MPICC is not set)
endif

CPPFLAGS = -DVERSION=$(RELEASE_VERSION)
CFLAGS = -flto $(MPI_COMPILE_FLAGS)
LDFLAGS = $(MPI_LINK_FLAGS)

msgrate: msgrate.o

msgrate.o: msgrate.c

clean:
	rm -f *.gcda msgrate.o msgrate *~

