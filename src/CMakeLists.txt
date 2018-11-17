cmake_minimum_required(VERSION 2.8)
project(FSC)

set(CMAKE_CXX_FLAGS "-std=c++11 -O3 -Wall")

set(SOURCE_FILES main.cpp
        Types.h
        Config.h
        CommandParser.h
        FSCCommand.h
        Graph.h
        Graph.cpp
        GraphReader.h
        GraphReader.cpp
        GraphDB.h
        GraphDB.cpp
        Utility.h
        Utility.cpp
        FSC.cpp
        FSC.h
        SequenceElement.h
        TreeNode.h)
add_executable(SubgraphQuery ${SOURCE_FILES})