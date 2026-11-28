execute_process(COMMAND "${CLIENT}" RESULT_VARIABLE status
    OUTPUT_VARIABLE result ERROR_VARIABLE detail TIMEOUT 30)
string(STRIP "${result}" result)
if(NOT status STREQUAL "0" OR NOT result STREQUAL "79228162514264337593543950336")
    message(FATAL_ERROR "Packaged exact-key client exited ${status}, output '${result}'\n${detail}")
endif()
