find_package(UnitTest++ REQUIRED)

add_executable(opinion_unit_tests main_tests.cpp)

target_link_libraries(opinion_unit_tests
                        PRIVATE
                        UnitTest++
                        opinion_model_lib
                        ${SOCP_LIBRARIES}
                      )

target_include_directories(opinion_unit_tests
                            PRIVATE
                            ${UTPP_INCLUDE_DIRS}
                            ${PROJECT_SOURCE_DIR}/include
                            ${SOCP_INCLUDE_DIR}
                          )
