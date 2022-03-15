#ifndef APERTURE_OPERATION_EXTRACTOR_H_
#define APERTURE_OPERATION_EXTRACTOR_H_


#include "independent_operation.hpp"

class Aperture_operation : public Independent_operation
{
  public:
    Aperture_operation(Lattice_element_slice const& slice)
      : Independent_operation("aperture")
    { }
};


class Finite_aperture_operation : public Aperture_operation
{
  public:
    Finite_aperture_operation(Lattice_element_slice const& slice)
      : Aperture_operation(slice) { }

  private:
    void apply_impl(Bunch & bunch, Logger & logger) const override { }

};

  std::unique_ptr<Independent_operation>
extract_aperture_operation(
    std::string const & aperture_type,
    Lattice_element_slice const & slice)
{
  if (aperture_type == "finite_aperture")
  {
    return std::make_unique<Finite_aperture_operation>(slice);
  }
  else
  {
    throw std::runtime_error("unknown aperture_type");
  }
}


#endif /* APERTURE_OPERATION_EXTRACTOR_H_ */

