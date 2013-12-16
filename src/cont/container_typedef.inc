         !----------------------------------------------------------------------
         ! Container class
         !----------------------------------------------------------------------
         TYPE, ABSTRACT :: ppm_t_container
             INTEGER :: nb = 0
             !!! number of different instances
             INTEGER :: iter_id = 0
             !!! index at which the iterator is currently pointing
             INTEGER :: min_id = 0
             !!! smallest non-null index in the array
             INTEGER :: max_id = 0
             !!! highest non-null index in the array
             INTEGER :: vec_size = 0
             !!! current allocated size of the array
         END TYPE ppm_t_container
