#ifndef AVT_GZ_BUFFERS_HH__
#define AVT_GZ_BUFFERS_HH__

#include "mesh_array.hh"
#include "loop_utils.hh"

namespace advect_1d {

    struct gz_buffer_t {
        private: 
        static std::size_t gz_size_ ;

        mesh_array_t<double> gz_data_ ;

    }

}

#endif 