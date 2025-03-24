



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
//// Helper function to calculate start and end indices:
////
array[] int calculate_start_indices_study( int n_studies,
                                           array[] int n_thr, 
                                           int n_tests) {
  
      array[n_tests] int start_index;
      
      start_index[1] = 1;
      for (t in 2:n_tests) {
        start_index[t] = start_index[t - 1] + (n_studies * n_thr[t - 1]);
      }
      
      return start_index;
  
}
////
//// Helper function to calculate start and end indices:
////
array[] int calculate_end_indices_study( int n_studies, 
                                         array[] int n_thr, 
                                         int n_tests,
                                         array[] int start_index) {
                                     
      array[n_tests] int end_index;
      
      for (t in 1:n_tests) {
        end_index[t] = start_index[t] + (n_thr[t] * n_studies) - 1;
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


  
  
  
  
  
/**
 * Gets a segment of values for a specific test and study from a vector
 *
 * @param original The vector containing all values
 * @param t The test index
 * @param s The study index
 * @param n_thr Array containing the number of thresholds for each test
 * @param n_studies Number of studies
 * @return Vector containing the segment for test t in study s
 */
vector get_test_t_study_s_segment_values( vector original, 
                                          int t, 
                                          int s, 
                                          array[] int n_thr, 
                                          int n_studies) {
                                            
        // Calculate the start index for this test-study combination
        int start_idx = 0;
        
        // Add all thresholds from previous tests across all studies
        for (prev_t in 1:(t-1)) {
          start_idx += n_thr[prev_t] * n_studies;
        }
        
        // Add thresholds for current test in previous studies
        start_idx += (s - 1) * n_thr[t];
        
        // Get the segment
        return segment(original, start_idx + 1, n_thr[t]);
  
}


/**
 * Updates a segment of a vector with new values
 *
 * @param flat_values_to_update The flat_values_to_update vector to be updated
 * @param new_values The new values to insert
 * @param t The test index
 * @param s The study index
 * @param n_thr Array containing the number of thresholds for each test
 * @param n_studies Number of studies
 * @return The updated vector
 */
vector update_test_t_study_s_segment_values( vector flat_values_to_update, 
                                             vector new_values, 
                                             int t, 
                                             int s, 
                                             array[] int n_thr, 
                                             int n_studies) {
                                                   
        vector[num_elements(flat_values_to_update)] result = flat_values_to_update;
        
        // Calculate the start index for this test-study combination
        int start_idx = 0;
        
        // Add all thresholds from previous tests across all studies
        for (prev_t in 1:(t-1)) {
          start_idx += n_thr[prev_t] * n_studies;
        }
        
        // Add thresholds for current test in previous studies
        start_idx += (s - 1) * n_thr[t];
        
        // Update the values
        for (i in 1:num_elements(new_values)) {
          result[start_idx + i] = new_values[i];
        }
        
        return result;
  
}
  








