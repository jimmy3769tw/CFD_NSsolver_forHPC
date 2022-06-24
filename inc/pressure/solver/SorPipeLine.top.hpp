
#pragma once


# include "controlPanel.hpp"

#define SOR_ERROR_MAXVAL
// #define SOR_ERROR_AVERAGE
// #define SOR_ERROR_SUM

#ifdef SOR_ERROR_AVERAGE
    #define SOR_ERROR_SUM
#endif
