find_package(GreatCMakeCookOff NO_MODULE PATHS ${cookoff_path} REQUIRED)
initialize_cookoff()
find_package(Julia)
include(CheckIsNaN)
include(TestCMake)

assert_recurse(-Dcookoff_path=${cookoff_path})
assert_recurse(FAIL -Dcookoff_path=${cookoff_path} -DDOFAILNOW=TRUE)
if(NORECURSE)
    if(NOT DOFAILNOW)
        find_package(Eigen)
    else()
        message(FATAL_ERROR "Should fail here")
    endif()
endif()
