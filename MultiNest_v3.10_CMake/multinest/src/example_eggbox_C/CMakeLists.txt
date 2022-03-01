cmake_minimum_required(VERSION 2.8) 

# In a normal application you would do:
# FIND_LIBRARY(MultiNest)
# INCLUDE_DIRECTORIES(${MultiNest_INCLUDE_PATH})

FIND_LIBRARY(m REQUIRED)

FILE(GLOB SOURCE *.c)
add_executable(eggboxC ${SOURCE})
target_link_libraries(eggboxC ${MultiNest_LIBRARIES} m)

# Set install directory
INSTALL(TARGETS eggboxC DESTINATION bin)
