#ifndef PARTICLE_TYPES_H
#define PARTICLE_TYPES_H

// We assume state[0–2] are positions and state[3–5] are momenta
struct ParticleState {
    double state[6];
    double t; // current simulation time
};

#endif // PARTICLE_TYPES_H
