         ! Write an attribute to a dataset
         INTERFACE write_attribute
            module procedure write_integer_attribute, &
                  write_double_attribute, &
                  write_logical_attribute, &
                  write_character_attribute, &
                  write_logical_array, &
                  write_character_array, &
                  write_string_attribute
         END INTERFACE write_attribute

         INTERFACE read_attribute
            MODULE PROCEDURE read_integer_attribute, &
                  read_double_attribute, &
                  read_string_attribute, &
                  read_logical_attribute, &
                  read_character_attribute, &
                  read_logical_array, &
                  read_character_array, &
                  read_logical_array_2d
         END INTERFACE read_attribute
