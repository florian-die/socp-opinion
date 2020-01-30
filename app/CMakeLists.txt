 add_executable(app_opinion testOpinion.cpp)

 target_link_libraries(app_opinion PRIVATE opinion_model_lib ${SOCP_LIBRARIES})
