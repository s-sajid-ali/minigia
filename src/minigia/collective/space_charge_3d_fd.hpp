#ifndef SPACE_CHARGE_3D_FD_H_
#define SPACE_CHARGE_3D_FD_H_

#include <minigia/bunch/core_diagnostics.hpp>
#include <minigia/simulation/collective_operator_options.hpp>
#include <minigia/simulation/operator.hpp>

#include "rectangular_grid_domain.hpp"

#include "space_charge_3d_fd_impl.hpp"

class Space_charge_3d_fd;

struct Space_charge_3d_fd_options
    : public CO_base_options<Space_charge_3d_fd_options, Space_charge_3d_fd> {
  std::array<int, 3> shape;
  bool domain_fixed;
  double n_sigma;
  double kick_scale;
  int comm_group_size;

  Space_charge_3d_fd_options(int gridx = 32, int gridy = 32, int gridz = 64)
      : shape{gridx, gridy, gridz}, domain_fixed(false), n_sigma(8.0),
        kick_scale(1.0), comm_group_size(1) {}

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
  const Space_charge_3d_fd_options options; /* Options to initialize sc-3d-fd */
  std::string bunch_sim_id;
  Rectangular_grid_domain domain;
  bool use_fixed_domain;
  bool allocated;

  GlobalCtx gctx;
  SubcommCtx sctx;
  LocalCtx lctx;

private:
  void apply_impl(Bunch_simulator &simulator, double time_step, Logger &logger);

  void get_local_charge_density(const Bunch &bunch);

  PetscErrorCode apply_bunch(Bunch &bunch, double time_step, Logger &logger);

  PetscErrorCode get_force();

  void apply_kick(Bunch &bunch, double time_step);

  PetscErrorCode allocate_sc3d_fd(const Bunch &bunch);

  PetscErrorCode destroy_sc3d_fd();

  PetscErrorCode update_domain(Bunch const &bunch);

public:
  Space_charge_3d_fd(Space_charge_3d_fd_options const &ops);

  void set_fixed_domain(std::array<double, 3> offset,
                        std::array<double, 3> size);

  ~Space_charge_3d_fd();
};

#endif /* SPACE_CHARGE_3D_FD_H_ */
