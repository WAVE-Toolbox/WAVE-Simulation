#pragma once

#define MASTER 0

/*!
 \def HOST_PRINT(comm,msg)
 Print message only on master.
 */
#define HOST_PRINT(comm,msg)     \
{                                   \
int myRank = comm->getRank();   \
if ( myRank == MASTER )         \
{                               \
std::cout << msg;           \
}                               \
}
