#ifndef STEPPER_H_
#define STEPPER_H_

#include <cereal/access.hpp>
#include <cereal/types/polymorphic.hpp>

#include "step.hpp"

class Stepper {
public:
  static const double fixed_step_tolerance;

public:
  virtual ~Stepper() = default;
  virtual std::unique_ptr<Stepper> clone() const = 0;

  std::vector<Step> apply(Lattice const &lattice) const {
    auto steps = apply_impl(lattice);
    create_operations(lattice, steps);
    return steps;
  }

private:
  void create_operations(Lattice const &lattice,
                         std::vector<Step> &steps) const {
    for (auto &step : steps)
      step.create_operations(lattice);
  }

  virtual std::vector<Step> apply_impl(Lattice const &lattice) const = 0;

  friend class cereal::access;

  template <class Archive> void serialize(Archive &ar) {}
};

// include the archive types before registering the derived class
#include <cereal/archives/json.hpp>

#endif /* STEPPER_H_ */
