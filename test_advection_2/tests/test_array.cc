#include<iostream>
#include<cstring>
#include<cstdlib>
#include<memory>
#include<vector>

template<typename T>
struct mesh_array_t {
    private:

    struct mesh_array_destructor {
        void operator()(void* buf){ free(buf) ; }
    };

    using avt_unique_ptr = std::unique_ptr<T, mesh_array_destructor> ;

    std::size_t size_ ;

    avt_unique_ptr data_ ;

    public: 

    mesh_array_t(): size_(0), data_(nullptr) {} ;

    mesh_array_t(std::size_t size__ ): size_(size__) {
        data_ = avt_unique_ptr( (T*) malloc(sizeof(T)*size_ )) ;
    }


    mesh_array_t(mesh_array_t<T>&& rhs_) : size_(rhs_.size_) {
        data_ = std::move(rhs_.data_) ;
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
        if( this == &rhs_ ){
            return *this ;
        }
        
        if( this->size_ != rhs_.size_ ) {
            data_.reset() ;
            size_ = rhs_.size_ ;
            data_ = avt_unique_ptr( (T*) malloc(sizeof(T)*size_ )) ;
        }

        std::copy( rhs_.data_.get(), rhs_.data_.get()+size_,this->data_.get()) ;
        return *this ;
    }

    std::size_t size() { return size_ ; }

    ~mesh_array_t() = default ;
};

template< bool have_rhs > 
struct field_t {
    private:
    std::size_t ntls_ ;
    std::size_t size_ ;

    std::vector<mesh_array_t<double>> data_ ;
    //mesh_array_t<double> data_ ;
    mesh_array_t<double> rhs_ ;

    public:

    field_t(
        std::size_t ntls__,
        std::size_t size__
    ):
    size_(size__),
    ntls_(ntls__)
    //data_(size_)
    {   
        
        data_.resize(ntls_) ;
        for( auto& x: data_) 
            x = mesh_array_t<double>(size_) ;
        
       /*
        for(int ii=0; ii<ntls_; ii++)
            data_.emplace_back(size_) ;
    */
        if( have_rhs )
            rhs_ = mesh_array_t<double>(size_) ;
    }

    
    mesh_array_t<double>& operator[](std::size_t tl_idx_ ){
        return data_[tl_idx_] ;
    }
    
    double& operator() ( std::size_t idx_ ){ return data_[0](idx_); }


};

int main() {

    field_t<false> field(2, 10) ;

    for(int ii=0; ii<10; ++ii) {
        field(ii) = 1. ;
        field[1](ii) = 2. ;
    }

    std::cout<< field[1](0) << std::endl ;

    mesh_array_t<double> a{2}, b{2} ;

    a(0) = a(1) = 1 ;
    a(1) = 1 ;
    b = a;

    std::cout << b(0) << std::endl ;
    std::cout << a(0) << std::endl ;


}