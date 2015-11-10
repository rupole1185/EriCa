project(MariSa)

file(GLOB MrsSOURCES *.c)

add_library( MariSa SHARED ${MrsSOURCES})

install(TARGETS MariSa DESTINATION lib)
install(DIRECTORY ./
        DESTINATION include/MariSa
        FILES_MATCHING
        PATTERN "*.h")
