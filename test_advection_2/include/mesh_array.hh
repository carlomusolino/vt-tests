#ifndef AVT_MESH_ARRAY_HH__
#define AVT_MESH_ARRAY_HH__

#include <vt/transport.h>
#include <memory>

#define AVT_ASSERT_SIZE_CHECK__

namespace advect_1d {

template<typename T>
struct mesh_array_t {
    private:

    struct mesh_array_destructor {
        void operator()(void* buf){  
            fmt::print("Destructor called\n");
            if( buf != nullptr  )
                vt::thePool()->dealloc(buf) ;
                //free(buf) ;
            }
    } ;

    using avt_unique_ptr = std::unique_ptr<T, mesh_array_destructor> ;
    std::size_t size_ ;

    avt_unique_ptr data_ ;

    public: 

    mesh_array_t(): size_(0), data_(nullptr)
    { } ;

    mesh_array_t(std::size_t size__ ): size_(size__) {
        fmt::print("Allocating buffer of size {}\n", size_) ;
        data_ = avt_unique_ptr( (T*) vt::thePool()->alloc(sizeof(T)*size_ )) ;
        //data_ = avt_unique_ptr( (T*) malloc(sizeof(T)*size_) ) ;
    }

    T& operator() (std::size_t idx) {
        #ifdef AVT_ASSERT_SIZE_CHECK__
        assert(idx < size_ ) ;
        #endif 
        return data_.get()[idx] ;
    }

    T const_access(std::size_t idx) const {
        #ifdef AVT_ASSERT_SIZE_CHECK__
        assert(idx < size_ ) ;
        #endif 
        return data_.get()[idx] ;
    }

    mesh_array_t<T>& operator=(mesh_array_t<T> const & rhs_ ){
        if( this == &rhs_ )
            return *this ;

        if ( this-> size_ != rhs_.size_ ) 
        {
            size_ = rhs_.size_ ;
            data_ .reset( (T*) vt::thePool()->alloc(size_ * sizeof(T) ) ) ;
        }
        // copy data
        std::copy( rhs_.data_.get(), 
            rhs_.data_.get()+size_, 
            this->data_.get() );
        
        return *this ;
    }

    mesh_array_t(mesh_array_t<T>&& rhs_ ): size_(rhs_.size_) {
        data_ = std::move(rhs_.data_) ;
    }
    /*
    mesh_array_t(const mesh_array_t<T>& rhs_ ): size_(rhs_.size_) {
        data_ = avt_unique_ptr( (T*) vt::thePool()->alloc(size_) ) ;
        std::copy( rhs_.data_.get(), 
            rhs_.data_.get()+size_, 
            this->data_.get() ) ;
    }*/

    void resize(const std::size_t & new_size__) {
        size_ = new_size__ ;
        data_.reset((T*) vt::thePool()->alloc(size_*sizeof(T))) ;
    }

    template< typename Serializer >
    void serialize(Serializer& s){
        s | size_ ;
        s | data_ ;
    }


    ~mesh_array_t() = default ;
};

}

#endif 