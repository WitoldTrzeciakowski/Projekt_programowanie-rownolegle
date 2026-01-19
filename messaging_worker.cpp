#include <sys/socket.h>
#include <netinet/in.h>
#include <unistd.h>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <algorithm>

#include "fitness.h"
#include "shared_types.h"

using namespace std;

double randDouble(double a, double b) {
    return a + (b - a) * (static_cast<double>(random()) / RAND_MAX);
}

double levyFlight() {
    const double u = randDouble(0, 1);
    const double v = randDouble(0, 1);
    return u / pow(fabs(v), 1.0 / 1.5);
}

void processLevyFlight(double* nests, double* fitnessVal, const int count, const int DIM,
                       const double LB, const double UB, const int fitness_id)
{
    const auto fitness = getFitness(fitness_id);
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

void processAbandonNests(double* nests, double* fitnessVal, const int count, const int DIM,
                         const double pa, const double LB, const double UB, const int fitness_id)
{
    const auto fitness = getFitness(fitness_id);
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

    listen(serverFd, 1);

    const int clientFd = accept(serverFd, nullptr, nullptr);

    while (true) {
        WorkPacket packet{};
        const ssize_t r = recv(clientFd, &packet, sizeof(packet), MSG_WAITALL);
        
        if (r <= 0 || packet.type == MSG_TERMINATE)
            break;

        const int count = packet.end - packet.start;
        constexpr int double_size = sizeof(double);
        const int nestBytes = count * packet.DIM * double_size;
        const int fitnessBytes = count * double_size;

        vector<double> nests(count * packet.DIM);
        vector<double> fitnessVal(count);

        recv(clientFd, nests.data(), nestBytes, MSG_WAITALL);
        recv(clientFd, fitnessVal.data(), fitnessBytes, MSG_WAITALL);

        if (packet.type == MSG_LEVY_FLIGHT) {
            processLevyFlight(nests.data(), fitnessVal.data(), count,
                              packet.DIM, packet.LB, packet.UB, packet.fitnessId);
        } else if (packet.type == MSG_ABANDON_NESTS) {
            processAbandonNests(nests.data(), fitnessVal.data(), count,
                                packet.DIM, packet.pa, packet.LB, packet.UB, packet.fitnessId);
        }

        send(clientFd, nests.data(), nestBytes, 0);
        send(clientFd, fitnessVal.data(), fitnessBytes, 0);
    }

    close(clientFd);
    close(serverFd);
    return 0;
}