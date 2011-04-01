#!/usr/bin/env python

""" This program prints out how to use a specific MPI system call """

usages = {
'mpi_init' : ['C  : MPI_Init(int argc, char* argv)',
              'F90: MPI_Init(ierror)',
              'Initializes the MPI environment'],
'mpi_finalize' : ['C  : MPI_Finalize()',
                  'F90: MPI_Finalize(ierror)',
                  'Finalizes MPI environment'],
'mpi_send' : ['C  : MPI_Send(*buffer, int count, MPI_Datatype, int destination, int tag, Communicator)',
              'F90: MPI_Send(buffer, int count, MPI_Datatype, int destination, int tag, Communicator, ierror)',
              'Sends data in a buffer from current thread to destination thread'],
'mpi_recv' : ['C  : MPI_Recv(*buffer, int count, MPI_Datatype, int source, int tag, Communicator), mpi_status',
              'F90: MPI_Recv(buffer, int count, MPI_Datatype, int source, int tag, Communicator, mpi_status, ierror)',
              'Receives data in a buffer from source thread'],
'mpi_comm_size' : ['C  : MPI_Comm_size(Communicator, int* size)',
                   'F90: MPI_Comm_size(Communicator, int size, ierror)',
                   'Returns the size of the communicator in "size"'],
'mpi_comm_rank' : ['C  : MPI_Comm_rank(Communicator, int* rank)',
                   'F90: MPI_Comm_rank(Communicator, int rank, ierror)',
                   'Returns the rank of the thread in "rank"'],
'mpi_bcast' : ['C  : MPI_Bcast(*message, int count, MPI_Datatype, int root, Communicator)',
               'F90: MPI_Bcast(message, int count, MPI_Datatype, int root, Communicator, ierror)',
               'Broadcasts data in "message" from root thread to every other thread in Communicator'],
'mpi_reduce' : ['C  : MPI_Reduce(*operand, *result, int count, MPI_Datatype, MPI_Operator, int root, Communicator)',
                'F90: MPI_Reduce(operand, result, int count, MPI_Datatype, MPI_Operator, int root, Communicator, ierror)',
                'Combines data in *operand into a *result on the root thread via the operator MPI_Operator.'],
'mpi_allreduce' : ['C  : MPI_Allreduce(*operand, *result, int count, MPI_Datatype, MPI_Operator, Communicator)',
                   'F90: MPI_Allreduce(operand, result, int count, MPI_Datatype, MPI_Operator, Communicator, ierror)',
                   'Combines data in *operand into a *result on all threads in Communicator'],
'mpi_gather' : ['C  : MPI_Gather(*sendbuf, int send_count, MPI_Datatype, *recvbuf, int recv_count, MPI_Datatype, int root, Communicator)',
                'F90: MPI_Gather(sendbuf, int send_count, MPI_Datatype, recvbuf, int recv_count, MPI_Datatype, int root, Communicator, ier)',
                'All threads dump sendbuf into recvbuf. Notes: send_count is how many values are sent from each thread; recv_count is how many\n' + \
                'values are received from each thread on root (NOT total received from all put together). However, recv_buffer must be allocated for all values.'],
'mpi_allgather' : ['C  : MPI_Gather(*sendbuf, int send_count, MPI_Datatype, *recvbuf, int recv_count, MPI_Datatype, Communicator)',
                   'F90: MPI_Gather(sendbuf, int send_count, MPI_Datatype, recvbuf, int recv_count, MPI_Datatype, Communicator, ier)',
                   'All threads dump sendbuf into recvbuf. Notes: send_count is how many values are sent from each thread; recv_count is how many\n' + \
                   'values are received from each thread on each thread (NOT total received from all put together). However, recv_buffer must be allocated for all values.'],
'mpi_scatter' : ['C  : MPI_Scatter(*sendbuf, int send_count, MPI_Datatype, *recvbuf, int recv_count, MPI_Datatype, Communicator)',
                 'F90: MPI_Scatter(sendbuf, int send_count, MPI_Datatype, recvbuf, int recv_count, MPI_Datatype, Communicator, ier)',
                 'Thread "root" splits sendbuf into numtasks different segments which it then distributes to each thread.'], 
'mpi_comm_split' : ['C  : MPI_Comm_split(Input Communicator, int color, int key, Output Communicator)',
                    'F90: MPI_Comm_split(Input Communicator, int color, int key, Output Communicator, ier)',
                    'Splits a communicator into several different communicators. "color" is the comm key you want to join\n' + \
                    '"key" is your rank in the input communicator'],
'mpi_barrier' : ['C  : MPI_Barrier(Communicator)', 'F90: MPI_Barrier(Communicator, ier)', 'Synchronizes threads by making them all wait']
}


import sys
from os import path

def printusage():
   """ Prints the usage and quits """
   print >> sys.stderr, 'Usage: %s <function>' % (path.split(sys.argv[0])[1])
   print >> sys.stderr, usages.keys()
   sys.exit(0)

if len(sys.argv) != 2 or sys.argv[1] == '--help' or sys.argv[1] == '-h':
   printusage()

# Now load a dictionary with usages and stuff

try:
   tmp = usages[sys.argv[1].lower()]
except:
   print >> sys.stderr, "%s not yet characterized!" % sys.argv[1]
   sys.exit(1)


print tmp[0]
print tmp[1]
print 'Description:\n%s' % tmp[2]
