#pragma once

// will only compile if AVX-512 (1st choice) or AVX2 (2nd choice) are available
// last vect_size elements
// Handle remainder
// If N < vect_size, handle everything with scalar operations
// if data in "long" format
// Row vector case
// Column vector case
// General matrix case
//// #pragma message "About to define AVX-512 mplementation of fn_process_Ref_double_AVX"
// fn_process_double_AVX_sub_function(x_Ref, 
//                                     test_simple_AVX512,  test_simple_double,  
//                                     test_simple_AVX512,  test_simple_double, skip_checks) ;
// fn_process_double_AVX_sub_function(x_Ref, test_simple_AVX512,  test_simple_double,   test_simple_AVX512, test_simple_double, skip_checks) ;
//// #pragma message "About to define AVX2 implementation of fn_process_Ref_double_AVX"
// fn_process_double_AVX_sub_function(x_Ref, 
//                                     test_simple_AVX2,  test_simple_double,  
//                                     test_simple_AVX2,  test_simple_double, skip_checks) ;
// fn_process_double_AVX_sub_function(x_Ref, test_simple_AVX2,  test_simple_double,   test_simple_AVX2, test_simple_double, skip_checks) ;
//// #pragma message "Defining dummy fn_process_Ref_double_AVX - since neither AVX2 nor AVX-512 are available"

#include <NicoSIMD/kernel_dispatch_and_vector_loops_of_the_BayesMVP_manual_gradient_models/fns_SIMD_and_wrappers/fn_wrappers_SIMD_AVX_general_for_MetaOrdDTA.hpp>























