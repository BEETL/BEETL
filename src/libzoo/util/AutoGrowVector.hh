#ifndef AUTO_GROW_VECTOR_HH
#define AUTO_GROW_VECTOR_HH

#include  <vector>


template <class T>
class AutoGrowVector : public std::vector<T>
{
public:
    T &operator[]( unsigned int n )
    {
        if ( this->size() <= n )
        {
            this->resize( n + 1 );
        }
        return ( ( std::vector<T> * )this )->operator[]( n );
    }

    const T &operator[]( unsigned int n ) const
    {
        if ( this->size() <= n )
        {
            this->resize( n + 1 );
        }
        return ( ( std::vector<T> * )this )->operator[]( n );
    }
};

#endif //AUTO_GROW_VECTOR_HH
