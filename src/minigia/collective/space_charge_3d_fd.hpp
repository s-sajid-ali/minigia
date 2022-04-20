#ifndef SPACE_CHARGE_3D_FD_H_
#define SPACE_CHARGE_3D_FD_H_

#include <minigia/bunch/core_diagnostics.hpp>
#include <minigia/simulation/collective_operator_options.hpp>
#include <minigia/simulation/operator.hpp>

#include "space_charge_3d_fd_impl.hpp"

class Space_charge_3d_fd;

struct Space_charge_3d_fd_options
    : public CO_base_options<Space_charge_3d_fd_options, Space_charge_3d_fd> {
  std::array<int, 3> shape;
  bool domain_fixed;
  double n_sigma;
  int comm_group_size;

  Space_charge_3d_fd_options(int gridx = 32, int gridy = 32, int gridz = 64)
      : shape{gridx, gridy, gridz}, domain_fixed(false), n_sigma(8.0),
        comm_group_size(1) {}

  template <class Archive> void serialize(Archive &ar) {
    ar(cereal::base_class<CO_base_options>(this));
    ar(shape);
    ar(n_sigma);
    ar(comm_group_size);
  }
};

CEREAL_REGISTER_TYPE(Space_charge_3d_fd_options);

/// Note: internal grid is stored in [z][y][x] order, but
/// grid shape expects [x][y][z] order.
class Space_charge_3d_fd : public Collective_operator {
private:
  const Space_charge_3d_fd_options options;
  std::string bunch_sim_id;
  bool use_fixed_domain;
  bool allocated;

  GlobalCtx gctx;
  SubcommCtx sctx;
  LocalCtx lctx;

private:
  void apply_impl(Bunch_simulator &simulator, double time_step, Logger &logger);

  void apply_bunch(Bunch &bunch, double time_step, Logger &logger);

  PetscErrorCode allocate_sc3d_fd(Space_charge_3d_fd_options const &ops);

  PetscErrorCode destroy_sc3d_fd();

  void update_domain(Bunch const &bunch);

public:
  Space_charge_3d_fd(Space_charge_3d_fd_options const &ops);
  ~Space_charge_3d_fd();
};

#endif /* SPACE_CHARGE_3D_FD_H_ */
