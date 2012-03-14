#define WRAP(a) a 
         TYPE, EXTENDS(WRAP(CONTAINER)_) ::  CONTAINER
                 CONTAINS
                     PROCEDURE       :: begin => __CONTAINER(begin)
                     PROCEDURE       :: next => __CONTAINER(next)
                     PROCEDURE       :: last => __CONTAINER(last)
                     PROCEDURE       :: prev => __CONTAINER(prev)
                     PROCEDURE       :: destroy => __CONTAINER(destroy)
                     PROCEDURE       :: exists => __CONTAINER(exists)
                     PROCEDURE       :: push => __CONTAINER(push)
                     PROCEDURE       :: remove => __CONTAINER(remove)
         END TYPE CONTAINER

#undef CONTAINER
#undef __CONTAINER
#undef VEC_TYPE
#undef WRAP
