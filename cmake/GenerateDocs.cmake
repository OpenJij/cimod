find_package(Doxygen
                 REQUIRED dot
                 OPTIONAL_COMPONENTS mscgen dia)
if (DOXYGEN_FOUND)
        include(GNUInstallDirs)
        find_program(DOXYGEN_EXECUTABLE doxygen REQUIRED)
        set(DOXYGEN_PROJECT_NAME "${PROJECT_NAME}")
        set(DOXYGEN_PROJECT_BRIEF "${PROJECT_DESCRIPTION}")
        set(DOXYGEN_EXTRACT_ALL YES)
        set(DOXYGEN_HAVE_DOT YES)
        set(DOXYGEN_DOT_MULTI_TARGETS YES)
        set(DOXYGEN_GENERATE_LATEX YES)
        set(DOXYGEN_RECURSIVE YES)
        set(DOXYGEN_GENERATE_HTML YES)
        set(DOXYGEN_CREATE_SUBDIRS YES)
        set(DOXYGEN_ALLOW_UNICODE_NAMES YES) 
        set(DOXYGEN_OUTPUT_LANGUAGE "Japanese-en")
        set(DOXYGEN_DOXYFILE_ENCODING "UTF-8")
        set(DOXYGEN_OUTPUT_TEXT_DIRECTION "ltr")
        set(DOXYGEN_OUTPUT_BRIEF_MEMBER_DESC YES)
        set(DOXYGEN_REPEAT_BRIEF YES)
        set(DOXYGEN_OUTPUT_
        set(DOXYGEN_OUTPUT_
        set(DOXYGEN_OUTPUT_
        set(DOXYGEN_OUTPUT_
        set(DOXYGEN_OUTPUT_
        set(DOXYGEN_OUTPUT_
        set(DOXYGEN_OPTIMIZE_OUTPUT_FOR_C YES)
        set(DOXYGEN_OPTIMIZE_OUTPUT_FOR_C YES)
        doxygen_add_docs(cxxcimod_header_only_docs
                         ${PROJECT_SOURCE_DIR}/src
                         ALL
                         COMMENT "Generate documentation with Doxygen")
        install(DIRECTORY ${PROJECT_BINARY_DIR}/docs/html
                DESTINATION ${CMAKE_INSTALL_DOCDIR}) 
    else() 
        message(SEND_ERROR "building documentation (-DENABLE_DOC=ON) is enabled, but doxygen not found")
    endif()
