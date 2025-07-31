


  
//// 
////  ---- Function to compute all pairwise differences between all of the thresholds of two tests
////  @param Se_or_Sp_t1 Vector of accuracies for all thresholds of test 1
////  @param Se_or_Sp_t2 Vector of accuracies for all thresholds of test 2
////  @param n_thr1 Number of thresholds for test 1
////  @param n_thr_t2 Number of thresholds for test 2
////  @return Vector of all pairwise differences
//// 
matrix compute_between_test_diffs( vector Se_or_Sp_t1, 
                                   vector Se_or_Sp_t2, 
                                   int n_thr_t1,
                                   int n_thr_t2) {
  
          matrix[n_thr_t1, n_thr_t2] pairwise_diffs;
          
          for (t1 in 1:n_thr_t1) {
            for (t2 in 1:n_thr_t2) {
              pairwise_diffs[t1, t2] = Se_or_Sp_t1[t1] - Se_or_Sp_t2[t2];
            }
          }
          
          return pairwise_diffs;
          
}



//// 
////  ---- Function to compute all pairwise ratios between all of the thresholds of two tests
////  @param Se_or_Sp_t1 Vector of accuracies for all thresholds of test 1
////  @param Se_or_Sp_t2 Vector of accuracies for all thresholds of test 2
////  @param n_thr1 Number of thresholds for test 1
////  @param n_thr_t2 Number of thresholds for test 2
////  @return Vector of all pairwise ratios
//// 
matrix compute_between_test_ratios( vector Se_or_Sp_t1, 
                                    vector Se_or_Sp_t2, 
                                    int n_thr_t1, 
                                    int n_thr_t2) {
  
        matrix[n_thr_t1, n_thr_t2] pairwise_ratios;
        
        for (t1 in 1:n_thr_t1) {
          for (t2 in 1:n_thr_t2) {
            pairwise_ratios[t1, t2] = Se_or_Sp_t1[t1] / Se_or_Sp_t2[t2];
          }
        }
        
        return pairwise_ratios;
  
}





































