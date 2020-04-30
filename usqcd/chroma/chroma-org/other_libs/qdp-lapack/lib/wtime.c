/* $Id: wtime.c,v 1.2 2008-04-22 03:53:05 kostas Exp $ */
 /*
 *   PRIMME PReconditioned Iterative MultiMethod Eigensolver
 *   Copyright (C) 2005  James R. McCombs,  Andreas Stathopoulos
 *
 *   This file is part of PRIMME.
 *
 *   PRIMME is free software; you can redistribute it and/or
 *   modify it under the terms of the GNU Lesser General Public
 *   License as published by the Free Software Foundation; either
 *   version 2.1 of the License, or (at your option) any later version.
 *
 *   PRIMME is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *   Lesser General Public License for more details.
 *
 *   You should have received a copy of the GNU Lesser General Public
 *   License along with this library; if not, write to the Free Software
 *   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */
#include <sys/time.h>
#include "sys/resource.h"
#ifdef RUSAGE_SELF
#else
#define   RUSAGE_SELF     0      /*needed in osx*/
#endif

double primme_wTimer(int zeroTimer) {
   static struct timeval tv;
   static struct timezone tz;
   static double StartingTime;
   
   if (zeroTimer) {
      gettimeofday(&tv, &tz); 
      StartingTime = ((double) tv.tv_sec) + ((double) tv.tv_usec )/(double) 1E6;
      return StartingTime;
   }
   else {
      gettimeofday(&tv, &tz); 
      return ((double) tv.tv_sec) + ((double) tv.tv_usec ) / (double) 1E6
	   - StartingTime;
   }
}

/* In the unlikely event that gettimeofday() is not available, but POSIX is, 
 * we can use the following alternative definition for primme_wTimer, 
 * after including time.h at the top.
 */
/*
#include <time.h>
double primme_wTimer(int zeroTimer) {
   static struct timespec ts;
   static double StartingTime;
   
   if (zeroTimer) {
      clock_gettime(CLOCK_REALTIME, &ts);
      StartingTime = ((double) ts.tv_sec) + ((double) ts.tv_nsec )/(double) 1E9;
      return StartingTime;
   }
   else {
      clock_gettime(CLOCK_REALTIME, &ts);
      return ((double) ts.tv_sec) + ((double) ts.tv_nsec ) / (double) 1E9;
	   - StartingTime;
   }
}
*/

/* 
 * Other timers that may be of use -------------------------------------------
 */

/* Simply return the microseconds time of day */
double primme_get_wtime() {
   static struct timeval tv;
   static struct timezone tz;

   gettimeofday(&tv, &tz); 
   return ((double) tv.tv_sec) + ((double) tv.tv_usec ) / (double) 1E6;
}

/* Return user/system times */
double primme_get_time(double *utime, double *stime) {
   struct rusage usage;
   static struct timeval utv,stv;

   getrusage(RUSAGE_SELF, &usage);
   utv = usage.ru_utime;
   stv = usage.ru_stime;

   *utime = ((double) utv.tv_sec) + ((double) utv.tv_usec ) / (double) 1E6;
   *stime = ((double) stv.tv_sec) + ((double) stv.tv_usec ) / (double) 1E6;

   return *utime + *stime;
}
