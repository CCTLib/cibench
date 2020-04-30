#ifndef SHIFT_TABLE_3D_PARSCALAR_H
#define SHIFT_TABLE_3D_PARSCALAR_H

#include <iostream>
#include <cstdlib>
#include <cache.h>
#include <qmp.h>

#include "shift_table_parscalar_types.h"

namespace CPlusPlusWilsonDslash {
  

  template< typename HalfSpinor >
  class ShiftTable3D {
  public: 
   
    ShiftTable3D(
	       const int* _subgrid_size,
	       HalfSpinor* chi1, 
	       HalfSpinor* chi2,
	       HalfSpinor* recv_bufs[2][3],
	       HalfSpinor* send_bufs[2][3],
	       void (*getSiteCoords)(int coord[], int node, int linearsite), 
	       int (*getLinearSiteIndex)(const int coord[]),
	       int (*getNodeNumber)(const int coord[]));

    ~ShiftTable3D() {
      free(xoffset_table);
      free(xsite_table);
    }

    inline
    int siteTable(int i) {
      return site_table[i];
    }

    HalfSpinor* halfspinorBufferOffset(HalfSpinorOffsetType type, int site, int mu) {
      //      std::cout << "type="<<type<<" site="<<site<<" mu="<<mu<<" index=" << (mu + 4*( site + subgrid_vol*(int)type)) << std::endl << std::flush;
	return offset_table[mu + 3*( site + subgrid_vol*(int)type) ];
    }

    inline int subgridVolCB() {
      return subgrid_vol_cb;
    }
  private:
    /* Tables */
    HalfSpinor** xoffset_table;        /* Unaligned */
    HalfSpinor** offset_table;         /* Aligned */
    
    int *xsite_table;         /* Unaligned */
    int *site_table;          /* Aligned */
        
    int tot_size[4];          /* Class scope members */
    int subgrid_size[4];
    int subgrid_cb_size[4];
    int subgrid_vol;
    int subgrid_vol_cb;         /* Useful numbers */
    const int Nd;                   /* No of Dimensions */
    

   // This is not needed as it can be done transitively:
  // ie lookup the QDP index and then lookup the coord with that 
    inline
    void mySiteCoords3D(int gcoords[], int node, int linearsite)
    {
      int mu;

      int tmp_coord[4];
      int cb3,cbb3;
      int* log_coords=QMP_get_logical_coordinates_from(node);
      int my_node = QMP_get_node_number();
      
      
      
      for(mu=0; mu < 4; mu++) { 
	gcoords[mu] = log_coords[mu]*subgrid_size[mu];
	
      }
      
      cb3=linearsite/subgrid_vol_cb;
      
      crtesn3d(linearsite % subgrid_vol_cb, subgrid_cb_size, tmp_coord);
      
      // Add on position within the node
      // NOTE: the cb for the x-coord is not yet determined
      gcoords[0] += 2*tmp_coord[0];
      for(mu=1; mu < 4; ++mu) {
	gcoords[mu] += tmp_coord[mu];
      }
      
      cbb3 = cb3;
      for(mu=1; mu < 3; ++mu) {
	cbb3 += gcoords[mu];
      }
      gcoords[0] += (cbb3 & 1);
    }
    
    int myLinearSiteIndex3D(const int gcoords[]) 
    {
      int mu;
      int subgrid_cb_coord[4];
      int cb3;
      const int Nd3 = 3;
      
      cb3=0;
      for(mu=0; mu < Nd3; ++mu) { 
	cb3 += gcoords[mu];
      }
      cb3 &=1;
      
      subgrid_cb_coord[0] = (gcoords[0]/2)% subgrid_cb_size[0];
      for(mu=1; mu < 4; mu++) { 
	subgrid_cb_coord[mu] = gcoords[mu] % subgrid_cb_size[mu];
      }

      return localSite4d(subgrid_cb_coord, subgrid_cb_size) + cb3*subgrid_vol_cb;
  }

    
    inline
    int localSite4d(const int coord[4], const int latt_size[4])
    {
      int order = 0;
      int mmu;
      
      // In the 4D Case: t+Lt(x + Lx(y + Ly*z)
      // essentially  starting from i = dim[Nd-2]
      //  order =  latt_size[i-1]*(coord[i])
      //   and need to wrap i-1 around to Nd-1 when it gets below 0
      for(mmu=2; mmu >= 0; --mmu) {
	int wrapmu = (mmu-1) % 4;
	if ( wrapmu < 0 ) wrapmu += 4;
	order = latt_size[wrapmu]*(coord[mmu] + order);
      }
      
      order += coord[ 3 ]; /* T is fastest running */
      
      return order;
    }
    
  
    inline
    void crtesn3d(int ipos, const int latt_size[], int coord[] )
    {
  
      int Ndim=3;
      int i, ix;
      
      /* Calculate the Cartesian coordinates of the VALUE of IPOS where the 
       * value is defined by
       *
       *     for i = 0 to NDIM-1  {
       *        X_i  <- mod( IPOS, L(i) )
       *        IPOS <- int( IPOS / L(i) )
       *     }
       *
       * NOTE: here the coord(i) and IPOS have their origin at 0. 
       */
      for(i = Ndim; i < Ndim+4; ++i) {
	ix=i%4;  /* This lets me start with the time direction and then wraparound */
	
	coord[ix] = ipos % latt_size[ix];
	ipos = ipos / latt_size[ix];
      }
      
    }


    inline 
    void offs(int temp[], const int coord[], int mu, int isign)
    {
      int i;
      
      for(i=0; i < 4; ++i) {
	temp[i] = coord[i];
      }
      /* translate address to neighbour */
      temp[mu] = (temp[mu] + isign + 2*tot_size[mu]) % tot_size[mu];

    } 

    inline
    int parity(const int coord[])
    {
      int m;
      int sum = 0;
      const int Nd3=3;

      for(m=0; m < Nd3; ++m)
	sum += coord[m];
      
      return sum & 1;
    }

  };



  // STUPID TEMPLATING MEANS I HAVE TO INLINE THIS...

#include  <cstddef>



 struct InvTab4 { 
   int cb;
   int linearcb;
 } ;
 

 template<typename HalfSpinor>
   ShiftTable3D<HalfSpinor>::ShiftTable3D(
				     const int* _subgrid_size,
				     HalfSpinor* chi1, 
				     HalfSpinor* chi2,
				     HalfSpinor* recv_bufs[2][3],
				     HalfSpinor* send_bufs[2][3],
				     void (*getSiteCoords)(int coord[], int node, int linearsite), 
				     int (*getLinearSiteIndex)(const int coord[]),
				     int (*getNodeNumber)(const int coord[])) : Nd(4)
  {


    /* Setup subgrid */
    const int* mach_size = QMP_get_logical_dimensions();
    int my_node = QMP_get_node_number();
    int bound[2][4][4];

    for(int mu=0; mu < Nd; mu++) { 
      subgrid_size[mu] = _subgrid_size[mu];
      subgrid_cb_size[mu] = _subgrid_size[mu];
      tot_size[mu] = mach_size[mu]*_subgrid_size[mu];
    }
    subgrid_cb_size[0] /= 2;

    subgrid_vol = subgrid_size[0];
    for(int mu=1; mu < 4; mu++) { 
      subgrid_vol *= subgrid_size[mu];
    }
    subgrid_vol_cb = subgrid_vol/2;


    /* Now I want to build the site table */
    /* I want it cache line aligned? */
    xsite_table = (int *)malloc(sizeof(int)*subgrid_vol+Cache::CacheLineSize);
    if(xsite_table == 0x0 ) { 
      QMP_error("Couldnt allocate site table");
      QMP_abort(1);
    }

    int pad = 0; 
    if( ((ptrdiff_t)xsite_table % Cache::CacheLineSize) != 0 ) { 
      pad = Cache::CacheLineSize - ((ptrdiff_t)xsite_table % Cache::CacheLineSize);
    }
    site_table = (int *)((unsigned char*)xsite_table+pad);

    /* I want an 'inverse site table'
       this is a one off, so I don't care so much about alignment 
     */
    InvTab4 *invtab = (InvTab4 *)malloc(sizeof(InvTab4)*subgrid_vol);
    if(invtab == 0x0 ) { 
      QMP_error("Couldnt allocate inv site table");
      QMP_abort(1);
    }


   /* Inversity of functions check:
       Check that myLinearSiteIndex3D is in fact the inverse
       of mySiteCoords3D, and that QDP_getSiteCoords is the
       inverse of QDP_linearSiteIndex()
    */
    for(int p=0; p < 2; p++) {
      for(int site=0; site < subgrid_vol_cb; site++) { 
	int gcoord[4];
	/* Linear site index */
	int my_index = site + subgrid_vol_cb*p;
	getSiteCoords(gcoord, my_node, my_index);
	int linear=getLinearSiteIndex(gcoord);

	if( linear != my_index ) { 
	  printf("P%d cb=%d site=%d : QDP_getSiteCoords not inverse of QDP_getLinearSiteIndex(): my_index=%d linear=%d\n", my_node, p,site, my_index,linear);
	}

	mySiteCoords3D(gcoord, my_node, my_index);
	linear=myLinearSiteIndex3D(gcoord);

	if( linear != my_index ) { 
	  printf("P%d cb=%d site=%d : mySiteCoords3D not inverse of myLinearSiteIndex3D(): my_index=%d linear=%d\n", my_node, p,site, my_index,linear);
	}
      }
    }

    const int* node_coord = QMP_get_logical_coordinates();

    for(int p=0; p < 2; p++) { 	      
      for(int z=0; z < subgrid_size[2]; z++) { 
	for(int y=0; y < subgrid_size[1]; y++) { 
	  for(int x=0; x < subgrid_size[0]/2; x++) { 
	    for(int t=0; t < subgrid_size[3]; t++) { 
	      
	      int coord[4];
	      coord[0] = 2*x + p;	      
	      coord[1] = y;
	      coord[2] = z;
	      coord[3] = t;

	      /* Make global */
	      for(int i=0; i < 4; i++) { 
		coord[i] += subgrid_size[i]*node_coord[i];
	      }

	      /* Both these indices serve as an index into something 
		 of lattice size. */

	      /* Index of coordinate -- NB this is not lexicographic
		 but takes into account checkerboarding in QDP++ */
	      int qdp_index = getLinearSiteIndex(coord);

	      /* Index of coordinate in my layout. -- NB this is not lexicographic
		 but takes into account my 3D checkerbaording */
	      int my_index = myLinearSiteIndex3D(coord);
	      

	      site_table[my_index] = qdp_index;

	      int cb3=parity(coord);
	      int linear = my_index%subgrid_vol_cb;

	      invtab[qdp_index].cb=cb3;
	      invtab[qdp_index].linearcb=linear;
	    }
	  }
	}
      }
    }


    
    /* Site table transitivity check: 
       for each site, convert to index in cb3d, convert to qdp index
       convert qdp_index to coordinate
       convert coordinate to back index in cb3d
       Check that your cb3d at the end is the same as you 
       started with */
    for(int p=0; p < 2; p++) { 
      for(int site=0; site < subgrid_vol_cb; site++) {
	int gcoord[4];

	/* My local index */
	int my_index = site + subgrid_vol_cb*p;
	
	/* Convert to QDP index */
	int qdp_index = site_table[ my_index ];
      
	/* Switch QDP index to coordinates */
	getSiteCoords(gcoord, my_node,qdp_index);
	
	/* Convert back to cb3d index */
	int linear = myLinearSiteIndex3D(gcoord);
	
	/* Check new cb,cbsite index matches the old cb index */
	if (linear != my_index) { 
	  printf("P%d The Circle is broken. My index=%d qdp_index=%d coords=%d,%d,%d,%d linear(=my_index?)=%d\n", my_node, my_index, qdp_index, gcoord[0],gcoord[1],gcoord[2],gcoord[3],linear);
	}
      }
    }


    for(int p=0; p < 2; p++) { 
      for(int site=0; site < subgrid_vol_cb; site++) {
      
	int gcoord[4];
	int gcoord2[4];

	/* My local index */
	int my_index = site + subgrid_vol_cb*p;
	mySiteCoords3D(gcoord, my_node, my_index);

	int qdp_index = site_table[ my_index ];
	getSiteCoords(gcoord2, my_node,qdp_index);
      
	for(int mu=0 ; mu < 4; mu++) { 
	  if( gcoord2[mu] != gcoord[mu] ) {
	    printf("P%d: my_index=%d qdp_index=%d mySiteCoords=(%d,%d,%d,%d) siteCoords=(%d,%d,%d,%d)\n", my_node, my_index, qdp_index, gcoord[0], gcoord[1], gcoord[2], gcoord[3], gcoord2[0], gcoord2[1], gcoord2[2], gcoord2[3]);
	    continue;
	  }
	}
      }
    }

    /* Allocate the shift table */
    /* The structure is as follows: There are 4 shift tables in order:
       
       [ Table 1 | Table 2 | Table 3 | Table 4 ]
       Table 1: decomp_scatter_index[mu][site]
       Table 2: decomp_hvv_scatter_index[mu][site]
       Table 3: recons_mvv_gather_index[mu][site]
       Table 4: recons_gather_index[mu][site]
       
    */
    int **shift_table;

    /* This 4 is for the 4 tables: Table 1-4*/
    if ((shift_table = (int **)malloc(4*sizeof(int*))) == 0 ) {
      QMP_error("init_wnxtsu3dslash: could not initialize shift_table");
    QMP_abort(1);
    
    }
  
    for(int i=0; i < 4; i++) { 
      /* This 4 is for the 4 comms dierctions: */
      if ((shift_table[i] = (int *)malloc(3*subgrid_vol*sizeof(int))) == 0) {
	QMP_error("init_wnxtsu3dslash: could not initialize shift_table");
	QMP_abort(1);
      }
    }

    /* Initialize the boundary counters */
    for(int cb=0; cb < 2; cb++) {
      for(int dir=0; dir < 3; dir++) {
	bound[cb][0][dir] = 0;	
	bound[cb][1][dir] = 0;	
	bound[cb][2][dir] = 0;	
	bound[cb][3][dir] = 0;	
      }
    }
    
    for(int cb=0; cb < 2; cb++) { 
      for(int site=0; site < subgrid_vol_cb; ++site) {
      
	int index = cb*subgrid_vol_cb + site;
      
	/* Fetch site from site table */
	int qdp_index = site_table[index];
      
	int coord[4];

	/* Get its coords */
	getSiteCoords(coord, my_node, qdp_index);
      
	/* Loop over directions building up shift tables */
	for(int dir=0; dir < 3; dir++) {
	
	  int fcoord[4], bcoord[4];
	  int fnode, bnode;
	  int blinear, flinear;
	
	  /* Backwards displacement*/
	  offs(bcoord, coord, dir, -1);
	  bnode   = getNodeNumber(bcoord);
	  blinear = getLinearSiteIndex(bcoord);

	  /* Forward displacement */
	  offs(fcoord, coord, dir, +1);
	  fnode   = getNodeNumber(fcoord);
	  flinear = getLinearSiteIndex(fcoord);
	  
	  /* Scatter:  decomp_{plus,minus} */
	  /* Operation: a^F(shift(x,type=0),dir) <- decomp(psi(x),dir) */ 
	  /* Send backwards - also called a receive from forward */
	  if (bnode != my_node) {      
	    /* Offnode */
	    /* Append to Tail 1, increase boundary count */
	    /* This is the correct code */
	    shift_table[DECOMP_SCATTER][dir+3*index] 
	      = subgrid_vol_cb + bound[1-cb][DECOMP_SCATTER][dir];
	    
	    bound[1-cb][DECOMP_SCATTER][dir]++;
	    
	  }
	  else {                                           
	    /* On node. Note the linear part of its (cb3, linear) bit,
	       using a reverse lookup */
	    shift_table[DECOMP_SCATTER][dir+3*index] = 
	      invtab[blinear].linearcb;
	  }
	
	
	  /* Scatter:  decomp_hvv_{plus,minus} */
	  /* Operation:  a^B(shift(x,type=1),dir) <- U^dag(x,dir)*decomp(psi(x),dir) */
	  /* Send forwards - also called a receive from backward */
	  if (fnode != my_node) {
	    /* Offnode */
	    /* Append to Tail 1, increase boundary count */
	    shift_table[DECOMP_HVV_SCATTER][dir+3*index]           
	      = subgrid_vol_cb + bound[1-cb][DECOMP_HVV_SCATTER][dir];
	    
	    bound[1-cb][DECOMP_HVV_SCATTER][dir]++;                  
	    
	  }
	  else {
	    /* On node. Note the linear part of its (cb3, linear) bit,
	       using a reverse lookup */
	    shift_table[DECOMP_HVV_SCATTER][dir+3*index]           /* Onnode */
	      = invtab[flinear].linearcb ;
	  }
	  
	  
	  /* Gather:  mvv_recons_{plus,minus} */
	  /* Operation:  chi(x) <-  \sum_dir U(x,dir)*a^F(shift(x,type=2),dir) */
	  /* Receive from forward */
	  if (fnode != my_node) {
	    /* Offnode */
	    /* Append to Tail 2, increase boundary count */
	    
	    shift_table[RECONS_MVV_GATHER][dir+3*index] =
	      2*subgrid_vol_cb + (bound[cb][RECONS_MVV_GATHER][dir]);
	    
	    bound[cb][RECONS_MVV_GATHER][dir]++;
	    
	  }
	  else {
	    /* On node. Note the linear part of its (cb3, linear) bit,
	       using a reverse lookup. Note this is a recons post shift,
	       so the linear coordinate to invert is mine rather than the neighbours */
	    shift_table[RECONS_MVV_GATHER][dir+3*index] =
	      invtab[qdp_index].linearcb ;
	  }
      
	  /* Gather:  recons_{plus,minus} */
	  /* Operation:  chi(x) +=  \sum_dir recons(a^B(shift(x,type=3),dir),dir) */
	  /* Receive from backward */
	  if (bnode != my_node) {
	    
	    shift_table[RECONS_GATHER][dir+3*index] = 
	      2*subgrid_vol_cb + bound[cb][RECONS_GATHER][dir];
	    
	    bound[cb][RECONS_GATHER][dir]++;
	    
	  }
	  else {
	    /* On node. Note the linear part of its (cb3, linear) bit,
	       using a reverse lookup. Note this is a recons post shift,
	       so the linear coordinate to invert is mine rather than the neighbours */
	  
	    shift_table[RECONS_GATHER][dir+3*index] = 
	      invtab[qdp_index].linearcb ;
	  }
	} 
      }
    }
    
    
    /* Sanity check - make sure the sending and receiving counters match */
    for(int cb=0; cb < 2; cb++) {
      for(int dir=0; dir < 3; dir++) {

	/* Sanity 1: Must have same number of boundary sites on each cb for 
	   a given operation */
	for(int i = 0; i < 4; i++) { 
	  if (bound[1-cb][i][dir] != bound[cb][i][dir]) {
	    
	    QMP_error("SSE Wilson dslash - make_shift_tables: type 0 diff. cb send/recv counts do not match: %d %d",
		      bound[1-cb][i][dir],bound[cb][i][dir]);
	    QMP_abort(1);
	  }
	}
      }
    }


      /* Now I want to make the offset table into the half spinor temporaries */
  /* The half spinor temporaries will look like this:
       
     dir=0 [ Body Half Spinors ][ Tail 1 Half Spinors ][ Tail 2 Half Spinors ]
     dir=1 [ Body Half Spinors ][ Tail 1 Half Spinors ][ Tail 2 Half Spinors ]
     ...
     
     And each of these blocks of half spinors will be sized to vol_cb
     sites (ie half volume only).  The shift_table() for a given site and
     direction indexes into one of these lines. So the offset table essentially
     delineates which line one picks, by adding an offset of 
     3*subgrid_vol_cb*dir 
     To the shift. The result from offset table, can be used directly as a
     pointer displacement on the temporaries. 
     
     Perhaps the best way to condsider this is to consider a value
     of shift_table[type][dir/site] that lands in the body. The
     shift table merely gives me a site index. But the data needs
     to be different for each direction for that site index. Hence 
     we need to replicate the body, for each dir. The 3xsubgrid_vol_cb
     is just there to take care of the buffers.

     Or another way to think of it is that there is a 'body element' index
     specified by the shift table lookup, and that dir is just the slowest
     varying index.
       
  */

  /* 4 dims, 4 types, rest of the magic is to align the thingie */
    xoffset_table = (HalfSpinor **)malloc(3*4*subgrid_vol*sizeof(HalfSpinor*)+Cache::CacheLineSize);
    if( xoffset_table == 0 ) {
      QMP_error("init_wnxtsu3dslash: could not initialize offset_table[i]");
      QMP_abort(1);
    }

    pad = 0; 
    if( ((ptrdiff_t)xoffset_table % Cache::CacheLineSize) != 0 ) { 
      pad = Cache::CacheLineSize - ((ptrdiff_t)xoffset_table % Cache::CacheLineSize);
    }    

    /* This is the bit what aligns straight from AMD Manual */
    offset_table = (HalfSpinor**)((unsigned char*)xoffset_table + pad);
    /* Walk through the shift_table and remap the offsets into actual
       pointers */

    /* DECOMP_SCATTER */
    int num=0;
    int offsite_found=0;
    int offset;
    for(int dir =0; dir < 3; dir++) { 

      /* Loop through all the sites. Remap the offsets either to local 
	 arrays or pointers */
      
      for(int site=0; site < subgrid_vol; site++) { 

	
	offset = shift_table[DECOMP_SCATTER][dir+3*site];
	 if( offset >= subgrid_vol_cb ) { 
	   /* Found an offsite guy. It's address must be to the send back buffer */
	   /* send to back index = recv from forward index = 0  */
	   offsite_found++;
	   
	   offset_table[ dir + 3*(site + subgrid_vol*DECOMP_SCATTER) ] =
	     send_bufs[0][num]+(offset - subgrid_vol_cb);
	 }
	 else { 
	   /* Guy is onsite: This is DECOMP_SCATTER so offset to chi1 */
	   offset_table[ dir + 3*(site + subgrid_vol*DECOMP_SCATTER) ] =
	     chi1+shift_table[DECOMP_SCATTER][dir+3*site]+subgrid_vol_cb*dir;
	 }
      }
      
      if( offsite_found > 0 ) { 
	/* If we found an offsite guy, next direction has to 
	   go into the next dir part of the send bufs */
	num++; 
      }
    }
    
    /* DECOMP_HVV_SCATTER */
    /* Restart num-s */
    num=0;
    for(int dir =0; dir <3; dir++) { 
      offsite_found=0;
      for(int site=0; site < subgrid_vol; site++) { 
	offset = shift_table[DECOMP_HVV_SCATTER][dir+3*site];
	if( offset >= subgrid_vol_cb ) { 
	  /* Found an offsite guy. It's address must be to the send forw buffer */
	  /* send to forward / receive from backward index = 1 */
	  offsite_found++;
	  
	  offset_table[ dir + 3*(site + subgrid_vol*DECOMP_HVV_SCATTER) ] =
	    send_bufs[1][num]+(offset - subgrid_vol_cb);
	}
	else { 
	  /* Guy is onsite. This is DECOMP_HVV_SCATTER so offset to chi2 */
	  offset_table[ dir + 3*(site + subgrid_vol*DECOMP_HVV_SCATTER) ] =
	    chi2+shift_table[DECOMP_HVV_SCATTER][dir+3*site ]+subgrid_vol_cb*dir;
	}
      }
      if( offsite_found > 0 ) { 
	num++; 
      }
    }

     /* RECONS_MVV_GATHER */
     num=0;
     for(int dir =0; dir <3; dir++) { 
       offsite_found=0;
       for(int site=0; site < subgrid_vol; site++) { 
	 offset = shift_table[RECONS_MVV_GATHER][dir+3*site];
	 if( offset >= 2*subgrid_vol_cb ) { 
	   /* Found an offsite guy. It's address must be to the recv from front buffer */
	   /* recv_from front index = send to back index = 0 */
	   offsite_found++;
	   offset_table[ dir + 3*(site + subgrid_vol*RECONS_MVV_GATHER) ] =
	     recv_bufs[0][num]+(offset - 2*subgrid_vol_cb);
	 }
	 else { 
	   /* Guy is onsite */
	   /* This is RECONS_MVV_GATHER so offset with respect to chi1 */
	   offset_table[ dir + 3*(site + subgrid_vol*RECONS_MVV_GATHER) ] =
	     chi1+shift_table[RECONS_MVV_GATHER][dir+3*site ]+subgrid_vol_cb*dir;
	 }
       }
       if( offsite_found > 0 ) { 
	 num++; 
       }
     }

     /* RECONS_GATHER */
     num=0;
     for(int dir=0; dir <3; dir++) { 
       offsite_found=0;
       for(int site=0; site < subgrid_vol; site++) { 
	 offset = shift_table[RECONS_GATHER][dir+3*site];
	 if( offset >= 2*subgrid_vol_cb ) { 
	   /* Found an offsite guy. It's address must be to the recv from back buffer */
	   /* receive from back = send to forward index =  1*/
	   offsite_found++;
	   offset_table[ dir + 3*(site + subgrid_vol*RECONS_GATHER) ] =
	     recv_bufs[1][num]+(offset - 2*subgrid_vol_cb);
	 }
	 else { 
	   /* Guy is onsite */
	   /* This is RECONS_GATHER so offset with respect to chi2 */
	   offset_table[ dir + 3*(site + subgrid_vol*RECONS_GATHER ) ] = 
	     chi2+shift_table[RECONS_GATHER][dir+3*site ]+subgrid_vol_cb*dir;
	 }
       }
       if( offsite_found > 0 ) { 
	 num++; 
       }
     }
     

     
     /* Free shift table - it is no longer needed. We deal solely with offsets */
     for(int i=0; i < 4; i++) { 
       free( (shift_table)[i] );
     }
     free( shift_table );
     
     free( invtab );
  }






  
} // Namespace

#endif
