cmake_minimum_required(VERSION 3.31)
project(Projekt_programowanie_rownolegle)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(rpc_server
        rpc_server.cpp
)


add_executable(messaging_worker
        messaging_worker.cpp
)


add_executable(cuckoo
        main.cpp
)
