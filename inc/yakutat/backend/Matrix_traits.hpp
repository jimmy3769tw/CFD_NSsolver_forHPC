#pragma once
#ifndef __MATRIX_TRAITS_HPP__
#define __MATRIX_TRAITS_HPP__

namespace yakutat::backend
{
    template<typename T>
    class MatrixTraits
    {
    public:
        using size_type = typename T::size_type;
        using value_type = typename T::value_type;
        using array_size_type = typename T::array_size_type;
        using array_value_type = typename T::array_value_type;
     
        
    };
}

#endif
