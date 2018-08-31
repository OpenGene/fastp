include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(fastp
    main.cpp
    adaptertrimmer.cpp
    basecorrector.cpp
    duplicate.cpp
    evaluator.cpp
    fastqreader.cpp
    filter.cpp
    filterresult.cpp
    htmlreporter.cpp
    jsonreporter.cpp
    nucleotidetree.cpp
    options.cpp
    overlapanalysis.cpp
    peprocessor.cpp
    polyx.cpp
    processor.cpp
    read.cpp
    seprocessor.cpp
    sequence.cpp
    stats.cpp
    threadconfig.cpp
    umiprocessor.cpp
    unittest.cpp
    writer.cpp
    writerthread.cpp)

target_link_libraries(fastp ${ZLIB_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

install(TARGETS fastp RUNTIME DESTINATION bin)
