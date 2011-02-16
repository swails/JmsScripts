#!/usr/bin/env python

""" This program prints out how to use a specific MPI system call """


import sys
from os import path

def printusage():
   """ Prints the usage and quits """
   print >> sys.stderr, 'Usage: %s <function>' % (path.split(sys.argv[0])[1])
   sys.exit(0)

if len(sys.argv) != 2 or sys.argv[1] == '--help' or sys.argv[1] == '-h':
   printusage()

# Now load a dictionary with usages and stuff

usages = {
'mpi_init' : ['C  : MPI_Init(int argc, char* argv)',
              'F90: MPI_Init(ierror)'],
'mpi_finalize' : ['C  : MPI_Finalize()',
                  'F90: MPI_Finalize(ierror)'],
'mpi_send' : ['C  : MPI_Send(*buffer, int count, MPI_Datatype, int destination, int tag, Communicator)',
              'F90: MPI_Send(buffer, int count, MPI_Datatype, int destination, int tag, Communicator, ierror)'],
'mpi_recv' : ['C  : MPI_Recv(*buffer, int count, MPI_Datatype, int source, int tag, Communicator)',
              'F90: MPI_Recv(buffer, int count, MPI_Datatype, int source, int tag, Communicator, ierror)'],
'mpi_comm_size' : ['C  : MPI_Comm_size(Communicator, int* size)',
                   'F90: MPI_Comm_size(Communicator, int size, ierror)'],
'mpi_comm_rank' : ['C  : MPI_Comm_rank(Communicator, int* rank)',
                   'F90: MPI_Comm_rank(Communicator, int rank, ierror)']
}

try:
   tmp = usages[sys.argv[1].lower()]
except:
   print >> sys.stderr, "%s not yet characterized!" % sys.argv[1]
   sys.exit(1)


print tmp[0]
print tmp[1]
