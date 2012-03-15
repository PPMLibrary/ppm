         TYPE,ABSTRACT,EXTENDS(ppm_t_container) ::  CONTAINER
             CLASS(VEC_TYPE),DIMENSION(:), POINTER :: vec => NULL()
             !!! array with pointers to meshes
             CLASS(VEC_TYPE),             POINTER :: iterator => NULL()
             CONTAINS
             PROCEDURE(__CONTAINER(begin)),  DEFERRED :: begin
             PROCEDURE(__CONTAINER(next)),   DEFERRED :: next 
             PROCEDURE(__CONTAINER(prev)),   DEFERRED :: prev 
             PROCEDURE(__CONTAINER(last)),   DEFERRED :: last 
             PROCEDURE(__CONTAINER(destroy)),DEFERRED :: destroy 
             PROCEDURE(__CONTAINER(exists)), DEFERRED :: exists 
             PROCEDURE(__CONTAINER(push)),   DEFERRED :: push 
             PROCEDURE(__CONTAINER(remove)), DEFERRED :: remove 
         END TYPE CONTAINER

#undef CONTAINER
#undef __CONTAINER
#undef VEC_TYPE
