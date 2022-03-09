#include "hdf5_writer.hpp"

namespace
{
    template <typename T>
        void storage_order_writer_helper(hid_t file_ptr, std::string const& name, T const& data)
        {
            int storage_order;

            if (data.storage_order() == boost::c_storage_order())
            {
                storage_order = Hdf5_writer<T>::c_storage_order;
            }
            else if (data.storage_order() == boost::fortran_storage_order())
            {
                storage_order = Hdf5_writer<T>::fortran_storage_order;
            }
            else
            {
                throw(std::runtime_error("Hdf5_writer: unsupported storage order"));
            }

            Hdf5_writer<int> storage_order_writer(file_ptr, name + "_storage_order");
            storage_order_writer.write(storage_order);
        }
}
