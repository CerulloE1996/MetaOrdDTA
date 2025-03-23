



////
//// Helper function to calculate start and end indices:
////
array[] int calculate_start_indices( array[] int n_thr, 
                                     int n_tests) {
  
      array[n_tests] int start_index;
      
      start_index[1] = 1;
      for (t in 2:n_tests) {
        start_index[t] = start_index[t - 1] + n_thr[t - 1];
      }
      
      return start_index;
  
}


////
//// Helper function to calculate start and end indices:
////
array[] int calculate_end_indices( array[] int n_thr, 
                                   int n_tests,
                                   array[] int start_index) {
                                     
      array[n_tests] int end_index;
      
      for (t in 1:n_tests) {
        end_index[t] = start_index[t] + n_thr[t] - 1;
      }
      
      return end_index;
  
}


////
//// Function to access elements from a flattened vector representing a ragged array:
////
vector get_test_values( vector flat_values, 
                        data array[] int start_index, 
                        data array[] int end_index, 
                        data int test_index) {
  
      int n_elements = end_index[test_index] - start_index[test_index] + 1;
      vector[n_elements] result;
      
      for (i in 1:n_elements) {
        result[i] = flat_values[start_index[test_index] + i - 1];
      }
      
      return result;
  
}


////
//// Function to set values in the flattened vector:
////
vector update_test_values( vector flat_values_to_update, 
                           vector new_values, 
                           data array[] int start_index, 
                           data array[] int end_index, 
                           data int test_index) {
  
      int n_elements = end_index[test_index] - start_index[test_index] + 1;
      vector[num_elements(flat_values_to_update)] result = flat_values_to_update;
      
      for (i in 1:n_elements) {
        result[start_index[test_index] + i - 1] = new_values[i];
      }
      
      return result;
      
}


  
  
  

