#include <sys/socket.h>
#include <netinet/in.h>
#include <unistd.h>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <algorithm>
#include <functional>
#include "shared_types.h"
#include "fitness.h"

using namespace std;

double randDouble(double a, double b) {
    return a + (b - a) * (static_cast<double>(random()) / RAND_MAX);
}

double levyFlight() {
    const double u = randDouble(0, 1);
    const double v = randDouble(0, 1);
    return u / pow(fabs(v), 1.0 / 1.5);
}

void executeLevyFlight(double* nests, double* fitnessVal, int count, int DIM,
                       double LB, double UB,
                       const function<double(const vector<double>&)>& fitness)
{
    vector<double> candidate(DIM);

    for (int i = 0; i < count; i++) {
        for (int d = 0; d < DIM; d++) {
            candidate[d] = nests[i * DIM + d] + levyFlight();
            candidate[d] = max(LB, min(UB, candidate[d]));
        }

        const double f_new = fitness(candidate);
        if (f_new < fitnessVal[i]) {
            for (int d = 0; d < DIM; d++)
                nests[i * DIM + d] = candidate[d];
            fitnessVal[i] = f_new;
        }
    }
}

void executeAbandonNests(double* nests, double* fitnessVal, int count, int DIM,
                         double pa, double LB, double UB,
                         const function<double(const vector<double>&)>& fitness)
{
    vector<double> tmp(DIM);

    for (int i = 0; i < count; i++) {
        if (randDouble(0, 1) < pa) {
            for (int d = 0; d < DIM; d++) {
                tmp[d] = randDouble(LB, UB);
                nests[i * DIM + d] = tmp[d];
            }
            fitnessVal[i] = fitness(tmp);
        }
    }
}

int main(int argc, char** argv) {
    if (argc < 2)
        return 1;

    const int port = atoi(argv[1]);
    srandom(getpid());

    const int serverFd = socket(AF_INET, SOCK_STREAM, 0);

    const int opt = 1;
    setsockopt(serverFd, SOL_SOCKET, SO_REUSEADDR, &opt, sizeof(opt));

    sockaddr_in addr{};
    addr.sin_family = AF_INET;
    addr.sin_addr.s_addr = INADDR_ANY;
    addr.sin_port = htons(port);

    const int result = bind(serverFd, reinterpret_cast<sockaddr*>(&addr), sizeof(addr));
    if (result != 0)
        return result;

    listen(serverFd, 8);

    while (true) {
        const int clientFd = accept(serverFd, nullptr, nullptr);

        RPCRequest req{};
        recv(clientFd, &req, sizeof(req), MSG_WAITALL);

        if (req.function == RPC_SHUTDOWN) {
            close(clientFd);
            break;
        }

        auto fitness = getFitness(req.fitnessId);
        constexpr int double_size = sizeof(double);

        const int nestBytes = req.nestCount * req.DIM * double_size;
        const int fitnessBytes = req.nestCount * double_size;

        vector<double> nests(req.nestCount * req.DIM);
        vector<double> fitnessVal(req.nestCount);

        recv(clientFd, nests.data(), nestBytes, MSG_WAITALL);
        recv(clientFd, fitnessVal.data(), fitnessBytes, MSG_WAITALL);

        if (req.function == RPC_LEVY_FLIGHT) {
            executeLevyFlight(nests.data(), fitnessVal.data(), req.nestCount,
                              req.DIM, req.LB, req.UB, fitness);
        } else if (req.function == RPC_ABANDON_NESTS) {
            executeAbandonNests(nests.data(), fitnessVal.data(), req.nestCount,
                                req.DIM, req.pa, req.LB, req.UB, fitness);
        }

        send(clientFd, nests.data(), nestBytes, 0);
        send(clientFd, fitnessVal.data(), fitnessBytes, 0);

        close(clientFd);
    }

    close(serverFd);
    return 0;
}