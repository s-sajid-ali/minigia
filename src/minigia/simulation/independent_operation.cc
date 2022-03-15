#include <minigia/libff/ff_element.hpp>

#include "independent_operation.hpp"

void
LibFF_operation::apply_impl(Bunch & bunch, Logger & logger) const
{
  for (auto const & slice : slices)
    FF_element::apply(slice, bunch);
}

