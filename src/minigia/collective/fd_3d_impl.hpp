#ifndef SPACE_CHARGE_3D_FD_IMPL_H_
#define SPACE_CHARGE_3D_FD_IMPL_H_

#include <petscdmda.h>
#include <petscvec.h>

/* Local (per MPI rank) context */
struct LocalCtx {

  Vec seqphi;  /*! local seq vector */
  Vec seqrho;  /*! local seq vector */
  DM da3d_seq; /*! sequential DMDA to manage grid and vecs */
};

#endif
