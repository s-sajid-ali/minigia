#include "resume.hpp"

Resume::Content::Content(Bunch_simulator * bunch_simulator_ptr,
    Stepper_sptr stepper_sptr) :
  bunch_sptr(bunch_simulator_ptr->get_bunch_sptr()), stepper_sptr(
      stepper_sptr), lattice_sptr(
        stepper_sptr->get_lattice_simulator().get_lattice_sptr())
{
}


Resume::Resume(std::string const& checkpoint_dir) :
  checkpoint_dir(checkpoint_dir), propagator()
{
  remove_serialization_directory();
  symlink_serialization_directory(checkpoint_dir);
  binary_load(
      propagator,
      get_combined_path(checkpoint_dir,
        Propagator::propagator_archive_name).c_str());
  unlink_serialization_directory();
}

  void
Resume::set_checkpoint_period(int period)
{
  propagator.set_checkpoint_period(period);
}

int
Resume::get_checkpoint_period() const
{
  return propagator.get_checkpoint_period();
}

  void
Resume::set_new_checkpoint_dir(std::string const& directory_name)
{
  propagator.set_checkpoint_dir(directory_name);
}

std::string const&
Resume::get_new_checkpoint_dir() const
{
  return propagator.get_checkpoint_dir();
}

  void
Resume::set_checkpoint_with_xml(bool with_xml)
{
  propagator.set_checkpoint_with_xml(with_xml);
}

bool
Resume::get_checkpoint_with_xml() const
{
  return propagator.get_checkpoint_with_xml();
}

  void
Resume::set_final_checkpoint(bool final_checkpoint)
{
  propagator.set_final_checkpoint(final_checkpoint);
}

bool
Resume::get_final_checkpoint() const
{
  return propagator.get_final_checkpoint();
}

  void
Resume::set_concurrent_io(int max)
{
  propagator.set_concurrent_io(max);
}

int
Resume::get_concurrent_io() const
{
  return propagator.get_concurrent_io();
}

  Resume::Content
Resume::get_content()
{
  Propagator::State state(propagator.get_resume_state(checkpoint_dir));
  Content content(state.bunch_simulator_ptr, propagator.get_stepper_sptr());
  return content;
}

  void
Resume::propagate(bool new_num_turns, int num_turns, bool new_max_turns, int max_turns, bool new_verbosity,
    int verbosity)
{
  propagator.resume(checkpoint_dir, new_num_turns, num_turns, new_max_turns, max_turns, new_verbosity,
      verbosity);
}

