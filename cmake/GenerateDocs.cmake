MESSAGE (STATUS "Build Documentation")
FIND_PACKAGE (Doxygen
              REQUIRED dot
              OPTIONAL_COMPONENTS mscgen dia)
IF (DOXYGEN_FOUND)
    INCLUDE (GNUInstallDirs)
    FIND_PROGRAM (DOXYGEN_EXECUTABLE doxygen REQUIRED)
    SET (DOXYGEN_PROJECT_NAME "${PROJECT_NAME}")
    SET (DOXYGEN_PROJECT_BRIEF "${PROJECT_DESCRIPTION}")
    
    
    SET (DOXYGEN_RECURSIVE YES)
    SET (DOXYGEN_GENERATE_HTML YES)
    
    #        set(DOXYGEN_CREATE_SUBDIRS YES)
    
    SET (DOXYGEN_ALLOW_UNICODE_NAMES YES)
    SET (DOXYGEN_OUTPUT_LANGUAGE "Japanese-en")
    SET (DOXYGEN_DOXYFILE_ENCODING "UTF-8")
    
    #        set(DOXYGEN_OUTPUT_TEXT_DIRECTION "ltr")
    
    SET (DOXYGEN_BRIEF_MEMBER_DESC YES)
    SET (DOXYGEN_REPEAT_BRIEF YES)
    #set(DOXYGEN_INLINE_INHERITED_MEMB YES)
    SET (DOXYGEN_INHERIT_DOCS YES)
    SET (DOXYGEN_DISTRIBUTE_GROUP_DOC YES)
    SET (DOXYGEN_GROUP_NESTED_COMPOUNDS YES)
    #set(DOXYGEN_INLINE_GROUPED_CLASSES YES)
    #set(DOXYGEN_INLINE_SIMPLE_STRUCTS YES)
    
    #        set(DOXYGEN_HTML_DYNAMIC_MENUS YES)
    
    SET (DOXYGEN_CHM_INDEX_ENCODING "UTF-8")
    
    #        set(DOXYGEN_HTML_DYNAMIC_SECTIONS YES)
    SET (DOXYGEN_HTML_TIMESTAMP YES)
    SET (DOXYGEN_ALWAYS_DETAILED_SEC YES)
    SET (DOXYGEN_JAVADOC_AUTOBRIEF YES)
    SET (DOXYGEN_FULL_PATH_NAMES NO)
    SET (DOXYGEN_SHORT_NAMES NO)
    SET (DOXYGEN_MULTILINE_CPP_IS_BRIEF NO)
    SET (DOXYGEN_INHERIT_DOCS YES)
    
    #        set(DOXYGEN_SEPARATE_MEMBER_PAGES YES)
    SET (DOXYGEN_TAB_SIZE 4)
    
    SET (DOXYGEN_MARKDOWN_SUPPORT YES)
    SET (DOXYGEN_AUTOLINK_SUPPORT YES)
    SET (DOXYGEN_BUILTIN_STL_SUPPORT NO)
    SET (DOXYGEN_CPP_CLI_SUPPORT YES)
    SET (DOXYGEN_IDL_PROPERTY_SUPPORT YES)
    SET (DOXYGEN_SUBGROUPING YES)
    SET (DOXYGEN_EXTRACT_ALL YES)
    SET (DOXYGEN_EXTRACT_PRIVATE YES)
    SET (DOXYGEN_EXTRACT_PACKAGE YES)
    SET (DOXYGEN_EXTRACT_PRIV_VIRTUAL YES)
    SET (DOXYGEN_EXTRACT_STATIC YES)
    SET (DOXYGEN_EXTRACT_LOCAL_CLASSES YES)
    SET (DOXYGEN_EXTRACT_LOCAL_METHODS YES)
    SET (DOXYGEN_EXTRACT_ANON_NSPACES YES)
    SET (DOXYGEN_SHOW_GROUPED_MEMB_INC YES)
    
    #        set(DOXYGEN_FORCE_LOCAL_INCLUDES YES)
    
    SET (DOXYGEN_INLINE_INFO YES)
    SET (DOXYGEN_SORT_MEMBER_DOCS YES)
    SET (DOXYGEN_SHOW_USED_FILES YES)
    SET (DOXYGEN_SHOW_FILES YES)
    SET (DOXYGEN_SHOW_NAMESPACES YES)
    SET (DOXYGEN_INPUT_ENCODING "UTF-8")
    SET (DOXYGEN_REFERENCES_LINK_SOURCE YES)
    SET (DOXYGEN_SOURCE_TOOLTIPS YES)
    #        set(DOXYGEN_ALPHABETICAL_INDEX YES)
    SET (DOXYGEN_GENERATE_TREEVIEW YES)
    SET (DOXYGEN_USE_MATHJAX YES)
    
    #        set(DOXYGEN_MACRO_EXPANSION YES)
    #        set(DOXYGEN_EXPAND_ONLY_PREDEF YES)
    #        set(DOXYGEN_PREDEFINED "extern=//")
    #        set(DOXYGEN_EXPAND_AS_DEFINED "extern")
    
    SET (DOXYGEN_CALL_GRAPH YES)
    SET (DOXYGEN_CALLER_GRAPH YES)
    SET (DOXYGEN_DOT_GRAPH_MAX_NODES 10000)
    SET (DOXYGEN_MAX_DOT_GRAPH_DEPTH 1000)
    SET (DOXYGEN_DIR_GRAPH_MAX_DEPTH 25)
    SET (DOXYGEN_HAVE_DOT YES)
    SET (DOXYGEN_DOT_MULTI_TARGETS YES)
    SET (DOXYGEN_CLASS_DIAGRAMS YES)
    SET (DOXYGEN_CLASS_GRAPH YES)
    SET (DOXYGEN_COLLABORATION_GRAPH YES)
    SET (DOXYGEN_DIRECTORY_GRAPH YES)
    SET (DOXYGEN_INCLUDE_GRAPH YES)
    SET (DOXYGEN_INCLUDED_BY_GRAPH YES)
    #        set(DOXYGEN_DOT_CLEANUP YES)
    
    #set(DOXYGEN_DOT_TRANSPARENT YES)
    SET (DOXYGEN_DOT_UML_DETAILS YES)
    SET (DOXYGEN_GRAPHICAL_HIERARCHY YES)
    SET (DOXYGEN_GROUP_GRAPHS YES)
    SET (DOXYGEN_INCLUDE_GRAPH YES)
    SET (DOXYGEN_INCLUDED_BY_GRAPH YES)
    #set(DOXYGEN_INTERACTIVE_SVG YES)
    SET (DOXYGEN_REFERENCED_BY_RELATION YES)
    SET (DOXYGEN_REFERENCES_RELATION YES)
    SET (DOXYGEN_GENERATE_XML YES)
    
    #set(DOXYGEN_UML_LIMIT_NUM_FIELDS 100)
    DOXYGEN_ADD_DOCS (cxxcimod_header_only_docs
                      ${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/cimod
                      ALL
                      COMMENT "Generate documentation with Doxygen")
    INSTALL (DIRECTORY ${PROJECT_BINARY_DIR}/html
             DESTINATION ${PROJECT_SOURCE_DIR}/docs)
ELSE ()
    MESSAGE (SEND_ERROR "building documentation (-DBUILD_DOCS=ON) is enabled, but doxygen not found")
ENDIF ()
