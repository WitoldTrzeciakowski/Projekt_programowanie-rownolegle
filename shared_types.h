#ifndef SHARED_TYPES_H
#define SHARED_TYPES_H

enum FitnessId {
    FITNESS_SPHERE = 1,
    FITNESS_ROSENBROCK = 2,
    FITNESS_RASTRIGIN = 3,
    FITNESS_ACKLEY = 4,
    FITNESS_SCHWEFEL = 5
};

enum MessageType { 
    MSG_LEVY_FLIGHT = 1, 
    MSG_ABANDON_NESTS = 2, 
    MSG_RESULT = 3,
    MSG_TERMINATE = 4 
};

enum RPCFunction { 
    RPC_LEVY_FLIGHT = 1, 
    RPC_ABANDON_NESTS = 2, 
    RPC_SHUTDOWN = 3 
};

struct WorkPacket {
    int type;
    int start;
    int end;
    int DIM;
    int N;
    double pa;
    double LB;
    double UB;
    int fitnessId;
};

struct RPCRequest {
    int function;
    int nestCount;
    int DIM;
    double pa;
    double LB;
    double UB;
    int fitnessId;
};

#endif //SHARED_TYPES_H