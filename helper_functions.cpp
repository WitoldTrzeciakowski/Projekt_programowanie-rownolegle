#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <thread>
#include <random>
#include <functional>

#include "fitness.h"
#include "shared_types.h"

using namespace std;

// ===================== helpers =====================

// Thread-local random number generator
thread_local mt19937 rng(random_device{}() + hash<thread::id>{}(this_thread::get_id()));
thread_local uniform_real_distribution<double> dist(0.0, 1.0);

double randDouble(const double a, const double b) {
    return a + (b - a) * (static_cast<double>(random()) / RAND_MAX);
}

double levyFlight() {
    const double u = randDouble(0, 1);
    const double v = randDouble(0, 1);
    return u / pow(fabs(v), 1.0 / 1.5);
}

double randDoubleParallel(const double a, const double b) {
    return a + (b - a) * dist(rng);
}

double levyFlightParallel() {
    const double u = dist(rng);
    const double v = dist(rng);
    return u / pow(fabs(v), 1.0 / 1.5);
}

vector<double> cuckooSearch(const int N, const int DIM, const int MAX_ITER,
                            const double pa, const double LB, const double UB,
                            const function<double(const vector<double>&)>& fitness)
{
    vector<vector<double>> nests(N, vector<double>(DIM));
    vector<double> fitnessVal(N);

    // Inicjalizacja
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < DIM; j++)
            nests[i][j] = randDouble(LB, UB);
        fitnessVal[i] = fitness(nests[i]);
    }

    for (int iter = 0; iter < MAX_ITER; iter++) {

        // --- Lévy flights ---
        for (int i = 0; i < N; i++) {
            vector<double> candidate = nests[i];

            for (int d = 0; d < DIM; d++) {
                candidate[d] += levyFlight();
                candidate[d] = max(LB, min(UB, candidate[d]));
            }

            const double f_new = fitness(candidate);
            if (f_new < fitnessVal[i]) {
                nests[i] = candidate;
                fitnessVal[i] = f_new;
            }
        }

        // --- Porzucanie gniazd ---
        for (int i = 0; i < N; i++) {
            if (static_cast<double>(random()) / RAND_MAX < pa) {
                for (int d = 0; d < DIM; d++)
                    nests[i][d] = randDouble(LB, UB);
                fitnessVal[i] = fitness(nests[i]);
            }
        }
    }

    // Wybór najlepszego rozwiązania
    const int64_t best = min_element(fitnessVal.begin(), fitnessVal.end()) - fitnessVal.begin();
    return nests[best];
}

// ===================== PARALLEL VERSION =====================

void levyFlightWorker(vector<vector<double>>& nests, 
                      vector<double>& fitnessVal,
                      const int start, const int end, const int DIM,
                      const double LB, const double UB,
                      const function<double(const vector<double>&)>& fitness)
{
    for (int i = start; i < end; i++) {
        vector<double> candidate = nests[i];

        for (int d = 0; d < DIM; d++) {
            candidate[d] += levyFlightParallel();
            candidate[d] = max(LB, min(UB, candidate[d]));
        }

        const double f_new = fitness(candidate);
        if (f_new < fitnessVal[i]) {
            nests[i] = candidate;
            fitnessVal[i] = f_new;
        }
    }
}

void abandonNestsWorker(vector<vector<double>>& nests,
                        vector<double>& fitnessVal,
                        const int start, const int end, const int DIM,
                        const double pa, const double LB, const double UB,
                        const function<double(const vector<double>&)>& fitness)
{
    for (int i = start; i < end; i++) {
        if (randDoubleParallel(0, 1) < pa) {
            for (int d = 0; d < DIM; d++)
                nests[i][d] = randDoubleParallel(LB, UB);
            fitnessVal[i] = fitness(nests[i]);
        }
    }
}

vector<double> cuckooSearchParallel(const int N, int DIM, const int MAX_ITER,
                                     double pa, double LB, double UB,
                                     function<double(const vector<double>&)> fitness,
                                     int num_threads = 0)
{
    if (num_threads <= 0)
        num_threads = static_cast<int>(thread::hardware_concurrency());
    
    vector<vector<double>> nests(N, vector<double>(DIM));
    vector<double> fitnessVal(N);

    // Inicjalizacja (parallel)
    vector<thread> threads;
    const int chunk = N / num_threads;
    
    for (int t = 0; t < num_threads; t++) {
        int start = t * chunk;
        int end = (t == num_threads - 1) ? N : (t + 1) * chunk;
        
        threads.emplace_back([&, start, end]() {
            for (int i = start; i < end; i++) {
                for (int j = 0; j < DIM; j++)
                    nests[i][j] = randDoubleParallel(LB, UB);
                fitnessVal[i] = fitness(nests[i]);
            }
        });
    }
    
    for (auto& t : threads)
        t.join();
    threads.clear();

    // Main loop
    for (int iter = 0; iter < MAX_ITER; iter++) {

        // --- Lévy flights (parallel) ---
        for (int t = 0; t < num_threads; t++) {
            int start = t * chunk;
            int end = (t == num_threads - 1) ? N : (t + 1) * chunk;
            
            threads.emplace_back(levyFlightWorker, ref(nests), ref(fitnessVal),
                               start, end, DIM, LB, UB, fitness);
        }
        
        for (auto& t : threads)
            t.join();
        threads.clear();

        // --- Porzucanie gniazd (parallel) ---
        for (int t = 0; t < num_threads; t++) {
            int start = t * chunk;
            int end = (t == num_threads - 1) ? N : (t + 1) * chunk;
            
            threads.emplace_back(abandonNestsWorker, ref(nests), ref(fitnessVal),
                               start, end, DIM, pa, LB, UB, fitness);
        }
        
        for (auto& t : threads)
            t.join();
        threads.clear();
    }

    // Wybór najlepszego rozwiązania
    const int64_t best = min_element(fitnessVal.begin(), fitnessVal.end()) - fitnessVal.begin();
    return nests[best];
}
// ===================== PARALLEL VERSION – PROCESSES =====================
// Linux / Unix only

#include <sys/mman.h>
#include <sys/wait.h>
#include <unistd.h>

// ---------- shared memory helpers ----------

double* allocSharedDouble(const size_t n) {
    return static_cast<double*>(mmap(
        nullptr,
        n * sizeof(double),
        PROT_READ | PROT_WRITE,
        MAP_SHARED | MAP_ANONYMOUS,
        -1,
        0
    ));
}

inline double& NEST(double* nests, const int i, const int d, const int DIM) {
    return nests[i * DIM + d];
}

// ---------- worker: Lévy flights ----------

void levyFlightProcess(double* nests, double* fitnessVal,
                       const int start, const int end, const int DIM,
                       const double LB, const double UB,
                       const function<double(const vector<double>&)>& fitness)
{
    vector<double> candidate(DIM);

    for (int i = start; i < end; i++) {
        for (int d = 0; d < DIM; d++) {
            candidate[d] = NEST(nests, i, d, DIM) + levyFlight();
            candidate[d] = max(LB, min(UB, candidate[d]));
        }

        const double f_new = fitness(candidate);
        if (f_new < fitnessVal[i]) {
            for (int d = 0; d < DIM; d++)
                NEST(nests, i, d, DIM) = candidate[d];
            fitnessVal[i] = f_new;
        }
    }
}

// ---------- worker: abandon nests ----------

void abandonNestsProcess(double* nests, double* fitnessVal,
                         const int start, const int end, const int DIM,
                         const double pa, const double LB, const double UB,
                         const function<double(const vector<double>&)>& fitness)
{
    vector<double> tmp(DIM);

    for (int i = start; i < end; i++) {
        if (randDouble(0, 1) < pa) {
            for (int d = 0; d < DIM; d++) {
                tmp[d] = randDouble(LB, UB);
                NEST(nests, i, d, DIM) = tmp[d];
            }
            fitnessVal[i] = fitness(tmp);
        }
    }
}

// ---------- MAIN: cuckoo search (process-based) ----------

vector<double> cuckooSearchProcess(const int N, const int DIM, const int MAX_ITER,
                                   const double pa, const double LB, const double UB,
                                   const function<double(const vector<double>&)>& fitness,
                                   int num_proc)
{
    if (num_proc <= 0)
        num_proc = static_cast<int>(sysconf(_SC_NPROCESSORS_ONLN));

    double* nests = allocSharedDouble(N * DIM);
    double* fitnessVal = allocSharedDouble(N);

    // ----- initialization -----
    for (int i = 0; i < N; i++) {
        vector<double> tmp(DIM);
        for (int d = 0; d < DIM; d++) {
            tmp[d] = randDouble(LB, UB);
            NEST(nests, i, d, DIM) = tmp[d];
        }
        fitnessVal[i] = fitness(tmp);
    }

    const int chunk = N / num_proc;

    for (int iter = 0; iter < MAX_ITER; iter++) {

        // ----- Lévy flights -----
        for (int p = 0; p < num_proc; p++) {
            if (fork() == 0) {
                const int start = p * chunk;
                const int end = (p == num_proc - 1) ? N : (p + 1) * chunk;
                levyFlightProcess(nests, fitnessVal,
                                  start, end, DIM, LB, UB, fitness);
                _exit(0);
            }
        }
        while (wait(nullptr) > 0) {}

        // ----- abandon nests -----
        for (int p = 0; p < num_proc; p++) {
            if (fork() == 0) {
                const int start = p * chunk;
                const int end = (p == num_proc - 1) ? N : (p + 1) * chunk;
                abandonNestsProcess(nests, fitnessVal,
                                    start, end, DIM, pa, LB, UB, fitness);
                _exit(0);
            }
        }
        while (wait(nullptr) > 0);
    }

    // ----- best solution -----
    int best = 0;
    for (int i = 1; i < N; i++)
        if (fitnessVal[i] < fitnessVal[best])
            best = i;

    vector<double> result(DIM);
    for (int d = 0; d < DIM; d++)
        result[d] = NEST(nests, best, d, DIM);

    munmap(nests, N * DIM * sizeof(double));
    munmap(fitnessVal, N * sizeof(double));

    return result;
}

// ===================== PARALLEL VERSION – MESSAGING =====================
// Linux / Unix only

#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>

void sendNests(const int sock, const double* nests, const double* fitness, const int start, const int end, const int DIM) {
    const int count = end - start;
    send(sock, &nests[start * DIM], count * DIM * sizeof(double), 0);
    send(sock, &fitness[start], count * sizeof(double), 0);
}

void recvNests(const int sock, double* nests, double* fitness, const int start, const int end, const int DIM) {
    const int count = end - start;
    recv(sock, &nests[start * DIM], count * DIM * sizeof(double), MSG_WAITALL);
    recv(sock, &fitness[start], count * sizeof(double), MSG_WAITALL);
}

vector<double> cuckooSearchMessaging(const int N, const int DIM, const int MAX_ITER,
                                      const double pa, const double LB, const double UB,
                                      const int fitness_id,
                                      int num_workers = 4,
                                      int basePort = 9000)
{
    auto fitness = getFitness(fitness_id);
    vector<double> nestsFlat(N * DIM);
    vector<double> fitnessVal(N);

    for (int i = 0; i < N; i++) {
        vector<double> tmp(DIM);
        for (int d = 0; d < DIM; d++) {
            tmp[d] = randDouble(LB, UB);
            nestsFlat[i * DIM + d] = tmp[d];
        }
        fitnessVal[i] = fitness(tmp);
    }

    vector<pid_t> workerPids;
    for (int w = 0; w < num_workers; w++) {
        pid_t pid = fork();
        if (pid == 0) {
            string portArg = to_string(basePort + w);
            execl("./messaging_worker", "messaging_worker", portArg.c_str(), nullptr);
            _exit(1);
        }
        workerPids.push_back(pid);
    }

    usleep(100000);

    vector<int> sockets(num_workers);
    for (int w = 0; w < num_workers; w++) {
        sockets[w] = socket(AF_INET, SOCK_STREAM, 0);

        sockaddr_in addr{};
        addr.sin_family = AF_INET;
        addr.sin_port = htons(basePort + w);
        inet_pton(AF_INET, "127.0.0.1", &addr.sin_addr);

        while (connect(sockets[w], reinterpret_cast<sockaddr*>(&addr), sizeof(addr)) < 0)
            usleep(10000);
    }

    int chunk = N / num_workers;

    for (int iter = 0; iter < MAX_ITER; iter++) {

        for (int w = 0; w < num_workers; w++) {
            int start = w * chunk;
            int end = (w == num_workers - 1) ? N : (w + 1) * chunk;

            WorkPacket packet{MSG_LEVY_FLIGHT, start, end, DIM, N, pa, LB, UB, fitness_id};
            send(sockets[w], &packet, sizeof(packet), 0);
            sendNests(sockets[w], nestsFlat.data(), fitnessVal.data(), start, end, DIM);
        }

        for (int w = 0; w < num_workers; w++) {
            int start = w * chunk;
            int end = (w == num_workers - 1) ? N : (w + 1) * chunk;
            recvNests(sockets[w], nestsFlat.data(), fitnessVal.data(), start, end, DIM);
        }

        for (int w = 0; w < num_workers; w++) {
            int start = w * chunk;
            int end = (w == num_workers - 1) ? N : (w + 1) * chunk;

            WorkPacket packet{MSG_ABANDON_NESTS, start, end, DIM, N, pa, LB, UB, fitness_id};
            send(sockets[w], &packet, sizeof(packet), 0);
            sendNests(sockets[w], nestsFlat.data(), fitnessVal.data(), start, end, DIM);
        }

        for (int w = 0; w < num_workers; w++) {
            int start = w * chunk;
            int end = (w == num_workers - 1) ? N : (w + 1) * chunk;
            recvNests(sockets[w], nestsFlat.data(), fitnessVal.data(), start, end, DIM);
        }
    }

    for (int w = 0; w < num_workers; w++) {
        WorkPacket packet{MSG_TERMINATE, 0, 0, 0, 0, 0, 0, 0, 0};
        send(sockets[w], &packet, sizeof(packet), 0);
        close(sockets[w]);
    }

    for (auto pid : workerPids)
        waitpid(pid, nullptr, 0);

    int best = 0;
    for (int i = 1; i < N; i++)
        if (fitnessVal[i] < fitnessVal[best])
            best = i;

    vector<double> result(DIM);
    for (int d = 0; d < DIM; d++)
        result[d] = nestsFlat[best * DIM + d];

    return result;
}

// ===================== PARALLEL VERSION – RPC =====================
// Linux / Unix only

int rpcConnect(const int port) {
    const int sock = socket(AF_INET, SOCK_STREAM, 0);

    sockaddr_in addr{};
    addr.sin_family = AF_INET;
    addr.sin_port = htons(port);
    inet_pton(AF_INET, "127.0.0.1", &addr.sin_addr);

    while (connect(sock, reinterpret_cast<sockaddr*>(&addr), sizeof(addr)) < 0)
        usleep(10000);

    return sock;
}

void rpcCall(const int port, const int function,
             double* nests, double* fitness,
             const int start, const int end, const int DIM,
             const double pa, const double LB, const double UB,
             const int fitness_id)
{
    const int sock = rpcConnect(port);

    const int count = end - start;
    const RPCRequest req{function, count, DIM, pa, LB, UB, fitness_id};

    send(sock, &req, sizeof(req), 0);
    send(sock, &nests[start * DIM], count * DIM * sizeof(double), 0);
    send(sock, &fitness[start], count * sizeof(double), 0);

    recv(sock, &nests[start * DIM], count * DIM * sizeof(double), MSG_WAITALL);
    recv(sock, &fitness[start], count * sizeof(double), MSG_WAITALL);

    close(sock);
}

void rpcShutdown(const int port) {
    const int sock = rpcConnect(port);
    constexpr RPCRequest req{RPC_SHUTDOWN, 0, 0, 0, 0, 0};
    send(sock, &req, sizeof(req), 0);
    close(sock);
}

vector<double> cuckooSearchRPC(const int N, const int DIM, const int MAX_ITER,
                                const double pa, const double LB, const double UB,
                                int fitness_id,
                                const int num_servers = 4,
                                const int basePort = 9100)
{
    vector<double> nestsFlat(N * DIM);
    vector<double> fitnessVal(N);
    const auto fitness = getFitness(fitness_id);

    for (int i = 0; i < N; i++) {
        vector<double> tmp(DIM);
        for (int d = 0; d < DIM; d++) {
            tmp[d] = randDouble(LB, UB);
            nestsFlat[i * DIM + d] = tmp[d];
        }
        fitnessVal[i] = fitness(tmp);
    }

    vector<pid_t> serverPids;
    for (int s = 0; s < num_servers; s++) {
        pid_t pid = fork();
        if (pid == 0) {
            const string portArg = to_string(basePort + s);
            execl("./rpc_server", "rpc_server", portArg.c_str(), nullptr);
            _exit(1);
        }
        serverPids.push_back(pid);
    }

    usleep(100000);

    const int chunk = N / num_servers;

    for (int iter = 0; iter < MAX_ITER; iter++) {

        vector<thread> threads;

        for (int s = 0; s < num_servers; s++) {
            int start = s * chunk;
            int end = (s == num_servers - 1) ? N : (s + 1) * chunk;

            threads.emplace_back(rpcCall, basePort + s, RPC_LEVY_FLIGHT,
                                 nestsFlat.data(), fitnessVal.data(),
                                 start, end, DIM, pa, LB, UB, fitness_id);
        }
        for (auto& t : threads) t.join();
        threads.clear();

        for (int s = 0; s < num_servers; s++) {
            int start = s * chunk;
            int end = (s == num_servers - 1) ? N : (s + 1) * chunk;

            threads.emplace_back(rpcCall, basePort + s, RPC_ABANDON_NESTS,
                                 nestsFlat.data(), fitnessVal.data(),
                                 start, end, DIM, pa, LB, UB, fitness_id);
        }
        for (auto& t : threads) t.join();
    }

    for (int s = 0; s < num_servers; s++)
        rpcShutdown(basePort + s);

    for (const auto pid : serverPids)
        waitpid(pid, nullptr, 0);

    int best = 0;
    for (int i = 1; i < N; i++)
        if (fitnessVal[i] < fitnessVal[best])
            best = i;

    vector<double> result(DIM);
    for (int d = 0; d < DIM; d++)
        result[d] = nestsFlat[best * DIM + d];

    return result;
}
