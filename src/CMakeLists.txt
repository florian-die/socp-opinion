set(OPINION_LIB_SOURCES opinion_model.cpp 
                        opinion_utils.cpp 
                        opinion_functions.cpp
                        )

set(OPINION_LIB_HEADERS ${PROJECT_SOURCE_DIR}/include/opinion_model.hpp
                        ${PROJECT_SOURCE_DIR}/include/opinion_params.hpp
                        ${PROJECT_SOURCE_DIR}/include/opinion_functions.hpp
                        )
 
add_library(opinion_model_lib STATIC ${OPINION_LIB_SOURCES} ${OPINION_LIB_HEADERS})                                      

target_link_libraries(opinion_model_lib PRIVATE ${SOCP_LIBRARIES})

target_include_directories(opinion_model_lib PRIVATE ${PROJECT_SOURCE_DIR}/include ${SOCP_INCLUDE_DIR})