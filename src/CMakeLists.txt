 add_library(opinion_model_lib STATIC opinion.cpp)

 target_link_libraries(opinion_model_lib PRIVATE ${SOCP_LIBRARIES})

 target_include_directories(opinion_model_lib PUBLIC ${PROJECT_SOURCE_DIR}/include)
