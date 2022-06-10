message(STATUS "Build Documentation")
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
        set(DOXYGEN_RECURSIVE YES)
        set(DOXYGEN_GENERATE_HTML YES)
#        set(DOXYGEN_CREATE_SUBDIRS YES)
        set(DOXYGEN_ALLOW_UNICODE_NAMES YES) 
        set(DOXYGEN_OUTPUT_LANGUAGE "Japanese-en")
        set(DOXYGEN_DOXYFILE_ENCODING "UTF-8")
        set(DOXYGEN_OUTPUT_TEXT_DIRECTION "ltr")
        set(DOXYGEN_BRIEF_MEMBER_DESC YES)
        set(DOXYGEN_REPEAT_BRIEF YES)
        
        set(DOXYGEN_ALWAYS_DETAILED_SEC YES)
        
        set(DOXYGEN_FULL_PATH_NAMES YES)
        set(DOXYGEN_SHORT_NAMES NO)
        set(DOXYGEN_MULTILINE_CPP_IS_BRIEF NO)
        set(DOXYGEN_INHERIT_DOCS YES)
#        set(DOXYGEN_SEPARATE_MEMBER_PAGES YES)
        set(DOXYGEN_TAB_SIZE 4)
        set(DOXYGEN_MARKDOWN_SUPPORT YES)
        set(DOXYGEN_AUTOLINK_SUPPORT YES)
        set(DOXYGEN_BUILTIN_STL_SUPPORT NO)
        set(DOXYGEN_CPP_CLI_SUPPORT NO)
        set(DOXYGEN_IDL_PROPERTY_SUPPORT NO)
        set(DOXYGEN_SUBGROUPING YES)
        set(DOXYGEN_EXTRACT_ALL YES)
        set(DOXYGEN_EXTRACT_PRIVATE YES)
        set(DOXYGEN_EXTRACT_PRIV_VIRTUAL YES)
        set(DOXYGEN_EXTRACT_STATIC YES)
        set(DOXYGEN_EXTRACT_LOCAL_CLASSES YES)
        set(DOXYGEN_EXTRACT_LOCAL_METHODS YES)
        set(DOXYGEN_EXTRACT_ANON_NSPACES YES)
        set(DOXYGEN_SHOW_GROUPED_MEMB_INC YES)
        set(DOXYGEN_FORCE_LOCAL_INCLUDES YES)
        set(DOXYGEN_INLINE_INFO YES)
        set(DOXYGEN_SORT_MEMBER_DOCS YES)
        set(DOXYGEN_SHOW_USED_FILES YES)
        set(DOXYGEN_SHOW_FILES YES)
        set(DOXYGEN_SHOW_NAMESPACES YES)
        set(DOXYGEN_INPUT_ENCODING "UTF-8")
        set(DOXYGEN_REFERENCES_LINK_SOURCE YES) 
        set(DOXYGEN_SOURCE_TOOLTIPS YES)
        set(DOXYGEN_ALPHABETICAL_INDEX YES)
        set(DOXYGEN_GENERATE_TREEVIEW YES)
        set(DOXYGEN_USE_MATHJAX YES)
        set(DOXYGEN_INCLUDE_GRAPH YES)
        set(DOXYGEN_INCLUDED_BY_GRAPH YES)
        doxygen_add_docs(cxxcimod_header_only_docs
                         ${PROJECT_SOURCE_DIR}/src
                         ALL
                         COMMENT "Generate documentation with Doxygen")
        install(DIRECTORY ${PROJECT_BINARY_DIR}/docs/html
                DESTINATION ${PROJECT_SOURCE_DIR}/docs) 
else() 
        message(SEND_ERROR "building documentation (-DBUILD_DOCS=ON) is enabled, but doxygen not found")
endif()
