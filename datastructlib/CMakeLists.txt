cmake_minimum_required(VERSION 3.10)
find_package(MPI REQUIRED)
# set the project name
project(datastructlib VERSION 1.1 DESCRIPTION "C Data Structures Library")

add_library(datastructlib SHARED
  arith.c
  assert.c
  atom.c
  except.c
  list.c
  queue.c
  stack.c
  table.c
  memchk.c
  )
