#include <vector>
#include <stdio.h>
#include "inc/track_gen.hpp"
#include "inc/trackUCN.hpp"
#include "inc/fields_nate.h"
#include "inc/lyap.hpp"
#include "inc/geometry.hpp"
#include <cmath>
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <string>
#include <getopt.h>
#include <cstring>

#include "setup.h"

//#include "inc/geometry.hpp"

extern "C" {
    #include "inc/xorshift.h"
}

void writeFixedRes(std::ofstream &binfile, fixedResult res) {
    const size_t buff_len = sizeof(unsigned int) + 9*sizeof(double) + 3*sizeof(int) + 2*sizeof(double) + sizeof(unsigned int);
    char buf[buff_len];
    if(!binfile.is_open()) {
        fprintf(stderr, "Error! file closed\n");
        return;
    }
    *((unsigned int *)(&buf[0])) = buff_len - 2*sizeof(unsigned int);
    *((double *)(&buf[0] + sizeof(unsigned int))) = res.energy;
    *((double *)(&buf[0] + sizeof(unsigned int) + 1*sizeof(double))) = res.theta;
    *((double *)(&buf[0] + sizeof(unsigned int) + 2*sizeof(double))) = res.t;
    *((double *)(&buf[0] + sizeof(unsigned int) + 3*sizeof(double))) = res.settlingT;
    *((double *)(&buf[0] + sizeof(unsigned int) + 4*sizeof(double))) = res.ePerp;
    *((double *)(&buf[0] + sizeof(unsigned int) + 5*sizeof(double))) = res.x;
    *((double *)(&buf[0] + sizeof(unsigned int) + 6*sizeof(double))) = res.y;
    *((double *)(&buf[0] + sizeof(unsigned int) + 7*sizeof(double))) = res.z;
    *((double *)(&buf[0] + sizeof(unsigned int) + 8*sizeof(double))) = res.zOff;
    *((int *)(&buf[0] + sizeof(unsigned int) + 9*sizeof(double))) = res.nHit;
    *((int *)(&buf[0] + sizeof(unsigned int) + 9*sizeof(double) + 1*sizeof(int))) = res.nHitHouseLow;
    *((int *)(&buf[0] + sizeof(unsigned int) + 9*sizeof(double) + 2*sizeof(int))) = res.nHitHouseHigh;
    *((double *)(&buf[0] + sizeof(unsigned int) + 9*sizeof(double) + 3*sizeof(int))) = res.eStart;
    *((double *)(&buf[0] + sizeof(unsigned int) + 9*sizeof(double) + 3*sizeof(int) + 1*sizeof(double))) = res.deathTime;
    *((unsigned int *)(&buf[0] + sizeof(unsigned int) + 9*sizeof(double) + 3*sizeof(int) + 2*sizeof(double))) = buff_len - 2*sizeof(unsigned int);
    binfile.write(buf, buff_len);
}

trace readTrace(const char *xfile, const char *yfile, const char *zfile) {
    trace t;
    t.x = NULL;
    t.y = NULL;
    t.z = NULL;
    t.num = -1;
    
    const size_t buff_len = 1*8;
    char* buf = new char[buff_len];
    
    std::ifstream binfileX(xfile, std::ios::in | std::ios::binary);
    std::ifstream binfileY(yfile, std::ios::in | std::ios::binary);
    std::ifstream binfileZ(zfile, std::ios::in | std::ios::binary);
    if(!binfileX.is_open() || !binfileY.is_open() || !binfileZ.is_open()) {
        fprintf(stderr, "Error! Could not open files!\n");
        return t;
    }
    
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;
    
    while(!binfileX.eof()) {
        binfileX.read(buf, buff_len);
        if(binfileX.eof()) { //Breaks on last read of file (i.e. when 0 bytes are read and EOF bit is set)
            break;
        }
        x.push_back(*(double *)&buf[0]);
    }
    binfileX.close();
    
    while(!binfileY.eof()) {
        binfileY.read(buf, buff_len);
        if(binfileY.eof()) { //Breaks on last read of file (i.e. when 0 bytes are read and EOF bit is set)
            break;
        }
        y.push_back(*(double *)&buf[0]);
    }
    binfileY.close();
    
    while(!binfileZ.eof()) {
        binfileZ.read(buf, buff_len);
        if(binfileZ.eof()) { //Breaks on last read of file (i.e. when 0 bytes are read and EOF bit is set)
            break;
        }
        z.push_back(*(double *)&buf[0]);
    }
    binfileZ.close();
    
    if(z.size() != x.size() || z.size() != y.size()) {
        fprintf(stderr, "Error! Sample length mismatch!\n");
        return t;
    }
    
    t.x = new double[x.size()];
    t.y = new double[y.size()];
    t.z = new double[z.size()];
    
    for(int i = 0; i < x.size(); i++) {
        t.x[i] = x[i];
        t.y[i] = y[i];
        t.z[i] = z[i];
    }
    
    t.num = x.size();
    
    return t;
}

// Kernel: Each thread integrates one trajectory.
__global__ void trackerKernel(ParticleState* d_states, fixedResult* d_results, int nTraj, double dt, trace* d_tr) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx >= nTraj) return;
    
    // Copy the 6-element state from global memory into a local array.
    double state[6];
    for (int i = 0; i < 6; i++) {
        state[i] = d_states[idx].state[i];
    }
    double t = d_states[idx].t;
    // Call the tracker device function.
    fixedResult res = fixedEffDaggerHitTime(state, dt, *d_tr);
    res.t += t;
    d_results[idx] = res;
}

int main(int argc, char** argv) {
    int ierr = MPI_Init(&argc, &argv);
    int nproc, rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    
    int c;
    char fName[1024];
    fName[0] = '\0';
    double dt = 0.0;
    double nTraj = 0;
    double PENindex = -1;
    
    if (PENTRACKSOURCE) {
        static struct option long_options[] = {
            {"file", required_argument, 0, 'f'},
            {"dt", required_argument, 0, 'd'},
            {"ntraj", required_argument, 0, 'n'},
            {"PENindex", required_argument, 0, 'p'},
            {0, 0, 0, 0},
        };
        while ((c = getopt_long(argc, argv, "f:d:n:p", long_options, NULL)) != -1) {
            switch(c) {
                case 'f':
                    strncpy(fName, optarg, 1024-1);
                    fName[1024-1] = '\0';
                    break;
                case 'd':
                    dt = atof(optarg);
                    break;
                case 'n':
                    nTraj = atof(optarg);
                    break;
                case 'p':
                    PENindex = atof(optarg);
                    break;
                default:
                    break;
            }
        }
        if(fName[0]=='\0' || dt==0 || nTraj==0 || PENindex==-1) {
            fprintf(stderr, "Usage: ./arrival_time_fixed_eff --file=fName --dt=timestep --ntraj=N --PENindex=startingIndex \n");
            exit(1);
        }
    } else {
        static struct option long_options[] = {
            {"file", required_argument, 0, 'f'},
            {"dt", required_argument, 0, 'd'},
            {"ntraj", required_argument, 0, 'n'},
            {0, 0, 0, 0},
        };
        while ((c = getopt_long(argc, argv, "f:d:n:", long_options, NULL)) != -1) {
            switch(c) {
                case 'f':
                    strncpy(fName, optarg, 1024-1);
                    fName[1024-1] = '\0';
                    break;
                case 'd':
                    dt = atof(optarg);
                    break;
                case 'n':
                    nTraj = atoi(optarg);
                    break;
                default:
                    break;
            }
        }
        if(fName[0]=='\0' || dt==0 || nTraj==0) {
            fprintf(stderr, "Usage: ./arrival_time_fixed_eff --file=fName --dt=timestep --ntraj=N\n");
            exit(1);
        }
    }
    
    // Read the trace files.
    trace tr = readTrace(XFNAME, YFNAME, ZFNAME);
    if(tr.x == NULL || tr.y == NULL || tr.z == NULL) {
        fprintf(stderr, "Bad trace\n");
        return 2;
    }
    printf("num trace bins: %d\n", tr.num);
    
    char fNameRank[1024];
    snprintf(fNameRank, 1024, "%s%d", fName, rank);
    std::ofstream binfile(fNameRank, std::ios::out | std::ios::binary);
    if (!binfile.is_open()) {
        fprintf(stderr, "Unable to open output file\n");
        return 1;
    }
    
    initxorshift(INITBLOCK);
    
    // Generate initial states for this MPI rank.
    std::vector<ParticleState> h_states;
    int localTraj = nTraj / nproc;
    h_states.reserve(localTraj);
    if (PENTRACKSOURCE) {
        std::ifstream myfile("neutrons_init.out");
        std::string line;
        if (myfile.is_open()) {
            getline(myfile, line); // Skip header
            for (int i = 0; i < PENindex; i++) {
                getline(myfile, line);
            }
            for (int i = 0; i < nTraj; i++){
                if (i >= rank*(nTraj/nproc) && i < (rank+1)*(nTraj/nproc)) {
                    if(getline(myfile, line)) {
                        std::vector<double> stateVec = readPENTrack(line, ' ');
                        ParticleState ps;
                        for (int j = 0; j < 6; j++) ps.state[j] = stateVec[j];
                        ps.t = 0.0;
                        h_states.push_back(ps);
                    } else {
                        fprintf(stderr, "Not enough neutrons in neutrons_init.out \n");
                        return 2;
                    }
                } else {
                    getline(myfile, line); // skip line
                }
            }
        } else {
            fprintf(stderr, "Unable to open neutrons_init.out\n");
            return 2;
        }
    } else {
        for (int i = 0; i < nTraj; i++) {
            if (i >= rank*(nTraj/nproc) && i < (rank+1)*(nTraj/nproc)) {
                std::vector<double> stateVec = TRACKGENERATOR(tr);
                ParticleState ps;
                for (int j = 0; j < 6; j++) ps.state[j] = stateVec[j];
                ps.t = 0.0;
                h_states.push_back(ps);
            } else {
                TRACKGENERATOR(tr);
            }
        }
    }
    
    // Allocate device memory.
    ParticleState* d_states;
    fixedResult* d_results;
    trace* d_tr;
    cudaMalloc(&d_states, localTraj * sizeof(ParticleState));
    cudaMalloc(&d_results, localTraj * sizeof(fixedResult));
    cudaMalloc(&d_tr, sizeof(trace));
    cudaMemcpy(d_tr, &tr, sizeof(trace), cudaMemcpyHostToDevice);
    cudaMemcpy(d_states, h_states.data(), localTraj * sizeof(ParticleState), cudaMemcpyHostToDevice);
    
    int threadsPerBlock = 256;
    int nBlocks = (localTraj + threadsPerBlock - 1) / threadsPerBlock;
    trackerKernel<<<nBlocks, threadsPerBlock>>>(d_states, d_results, localTraj, dt, d_tr);
    cudaDeviceSynchronize();
    
    std::vector<fixedResult> h_results(localTraj);
    cudaMemcpy(h_results.data(), d_results, localTraj * sizeof(fixedResult), cudaMemcpyDeviceToHost);
    
    // Write results using your WRITER (assumed to be defined elsewhere).
    for (auto &res : h_results) {
        WRITER(binfile, res);
    }
    
    binfile.close();
    cudaFree(d_states);
    cudaFree(d_results);
    cudaFree(d_tr);
    MPI_Finalize();
    return 0;
}
