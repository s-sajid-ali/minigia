
#include <minigia/bunch/period.hpp>
#include <minigia/foundation/physical_constants.hpp>

#include "step.hpp"

namespace {
void apply_longitudinal_boundary(Bunch &bunch) {
  // Bunch longitudinal boundary condition
  auto lb = bunch.get_longitudinal_boundary();

  switch (lb.first) {
  case LongitudinalBoundary::periodic:
    apply_longitudinal_periodicity(bunch, lb.second);
    break;

  case LongitudinalBoundary::aperture:
    apply_zcut(bunch, lb.second);
    break;

  case LongitudinalBoundary::bucket_barrier:
    apply_longitudinal_bucket_barrier(bunch, lb.second);
    break;

  case LongitudinalBoundary::open:
  default:
    break;
  }
}
} // namespace

Step::Step(double length) : operators(), length(length), step_betas() {}

void Step::create_operations(Lattice const &lattice) {
  for (auto &op : operators) {
    op->create_operations(lattice);
  }
}

void Step::apply(Bunch_simulator &simulator, Logger &logger) const {
  if (simulator[0].get_bunch_array_size() == 0) {
    throw std::runtime_error(
        "Step::apply() unable to proceed. no bunch in the simulator");
  }

  // time [s] in accelerator frame
  double ref_beta = simulator[0][0].get_reference_particle().get_beta();
  double time = length / (ref_beta * pconstants::c);

  for (auto const &op : operators) {
    double t0 = MPI_Wtime();

    logger(LoggerV::INFO_OPR) << "\n  Operator start:\n";

    // operator apply
    op->apply(simulator, time, logger);

    double t1 = MPI_Wtime();

    logger(LoggerV::INFO_OPR)
        << "  Operator finish: operator: name = " << op->get_name()
        << ", type = " << op->get_type() << ", time = " << std::fixed
        << std::setprecision(3) << t1 - t0 << "s"
        << "\n";

    // per operator diagnostics action
    simulator.diag_action_operator(*op);

    // longitudinal conditions
    for (auto &train : simulator.get_trains())
      for (auto &bunch : train.get_bunches())
        apply_longitudinal_boundary(bunch);
  }
}

#if 0
  void
Step::set_betas(double betax, double betay)
{
  this->step_betas.push_back(betax);
  this->step_betas.push_back(betay);
}

  std::vector<double >
Step::get_betas()
{
  return step_betas;
}

void
Step::print(int index) const
{
#if 0
  std::cout << "step " << index << ":\n";
  for (Operators::const_iterator it = operators.begin(); it
      != operators.end(); ++it) {
    (*it)->print();
  }
#endif
}

#endif
