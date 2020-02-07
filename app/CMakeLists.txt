 add_executable(app_opinion main_app.cpp)

 target_link_libraries(app_opinion PRIVATE opinion_model_lib ${SOCP_LIBRARIES})

 target_include_directories(app_opinion PRIVATE ${PROJECT_SOURCE_DIR}/include ${SOCP_INCLUDE_DIR})
