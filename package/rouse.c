#include "mex.h"
#include "math.h"
#include "stdlib.h"
#include "string.h"

//https://stackoverflow.com/questions/2325472/generate-random-numbers-following-a-normal-distribution-in-c-c
// double sampleNormal() {
//     double u = ((double) rand() / (RAND_MAX)) * 2 - 1;
//     double v = ((double) rand() / (RAND_MAX)) * 2 - 1;
//     double r = u * u + v * v;
//     if (r == 0 || r > 1) return sampleNormal();
//     double c = sqrt(-2 * log(r) / r);
//     return u * c;
// }

/* Function to generate a single random number with a normal distribution */
double generate_normal_random(double mean, double stddev) {
    static const double TWO_PI = 6.28318530717958647692;
    static int hasSpare = 0;
    static double spare;
    
    if (hasSpare) {
        hasSpare = 0;
        return mean + stddev * spare;
    }
    
    hasSpare = 1;
    
    double u1, u2;
    do {
        u1 = rand() / ((double)RAND_MAX);
        u2 = rand() / ((double)RAND_MAX);
    } while (u1 <= 1e-10);  // Avoid log(0)

    double r = sqrt(-2.0 * log(u1));
    double theta = TWO_PI * u2;

    spare = r * sin(theta);
    return mean + stddev * (r * cos(theta));
}

// MEX function definition
void rouse_sim(double Xs, int Nb, double zeta, double ks, double kT, double h, double fTime, int frameLen, int numFrames, double *sol) {

    // Input parameters equivalent to MATLAB
    int steps = (int)(fTime / h);  // Total number of steps
    
    // Solution output: only the frames we want..
    double *rmatCur = (double*)malloc(Nb * sizeof(double));  // Bead positions
    double *rmatNext = (double*)malloc(Nb * sizeof(double)); // Next bead positions
    double *rMat = (double*)malloc(Nb * Nb * sizeof(double)); // Coefficient matrix
    double *constArr = (double*)malloc(Nb * sizeof(double)); // Constants array
    double *curRnd = (double*)malloc(Nb * sizeof(double));   // Random noise vector
    double *f1 = (double*)malloc(Nb * sizeof(double));        // Force array f1
    double *f2 = (double*)malloc(Nb * sizeof(double));        // Force array f2
    double *initialGuess = (double*)malloc(Nb * sizeof(double)); // Initial guess for the next step

    // Initialize bead positions (rmatCur)
    for (int i = 0; i < Nb; i++) {
        rmatCur[i] = i * Xs;
        rmatNext[i] = 0;
        constArr[i] = 0;
        curRnd[i] = 0;
        f1[i] = 0;
        f2[i] = 0;
        initialGuess[i] = 0;
    }

//     printf("Initial rmatCur:\n");
//     for (int i = 0; i < Nb; i++) {
//             printf("%.8f ", rmatCur[i]);
//     }
    
    // Boundary conditions (constant positions)
    constArr[0] = -Xs;
    constArr[Nb - 1] = Xs;

    // debugging: print constArr
        // Debugging: Print the matrix rMat to verify its values
//     printf("Initial constArr:\n");
//     for (int i = 0; i < Nb; i++) {
//             printf("%.8f ", constArr[i]);
//     }
    
    // Initialize rMat as a tridiagonal matrix (NxN)
    for (int i = 0; i < Nb; i++) {
        for (int j = 0; j < Nb; j++) {
            if (i == j) {
                rMat[i * Nb + j] = -2; // Diagonal elements are -2
            } else if (i == j - 1 || i == j + 1) {
                rMat[i * Nb + j] = 1;  // Off-diagonal elements are 1
            } else {
                rMat[i * Nb + j] = 0;  // All other elements are 0
            }
        }
    }
    
    // Boundary conditions for the matrix (first and last rows)
    rMat[0 * Nb + 0] = -1;
    rMat[0 * Nb + 1] = 1;
    rMat[(Nb-1) * Nb + (Nb-2)] = 1;
    rMat[(Nb-1) * Nb + (Nb-1)] = -1;

    // Debugging: Print the matrix rMat to verify its values
//     printf("Initial rMat:\n");
//     for (int i = 0; i < Nb; i++) {
//         for (int j = 0; j < Nb; j++) {
//             printf("%.2f ", rMat[i * Nb + j]);
//         }
//         printf("\n");
//     }

    // Calculate the total number of time steps
    steps = (int)(fTime / h);

    // Compute the interval at which frames are stored
    int frameInterval = steps / numFrames;

    // Initialize solFrames to store the indices of the steps corresponding to frames
    int *solFrames = (int*)malloc(numFrames * sizeof(int));
    for (int i = 0; i < numFrames; i++) {
        solFrames[i] = (i + 1) * frameInterval;
    }

//     // Debugging: Print the solFrames to verify the indices
//     for (int i = 0; i < numFrames; i++) {
//         printf("Frame %d corresponds to step %d\n", i + 1, solFrames[i]);
//     }

    for (int i = 0; i < Nb; i++) {
        sol[i] = i * Xs;  // Equidistant positions at first frame
    }

    // Main simulation loop
    int k = 0;
    for (int t = 1; t <= steps; t++) {
        
        // Generate random noise (curRnd)
        for (int i = 0; i < Nb; i++) {
            curRnd[i] = sqrt(2 * kT / (zeta * h))* (generate_normal_random(0,1));  // Gaussian noise approximation
        }

//         if (k < numFrames && t >= solFrames[k]) {
//         // Store updated positions at subsequent frames
//         for (int i = 0; i < Nb; i++) {
//         sol[(k+1) * Nb + i] = curRnd[i];
//         }
//         k++;
//         }
//         for (int i = 0; i < Nb; i++) {
//             printf("%.8f ", curRnd[i]);
//         }

        f1[0] = (ks / zeta) * (rMat[0] * rmatCur[0]+rMat[1] * rmatCur[1]  + constArr[0]) + curRnd[0];
        // First step: Euler method (f1)
        for (int i = 1; i < Nb-1; i++) {
            f1[i] = (ks / zeta) * (rMat[i * Nb + i-1] * rmatCur[i-1] +rMat[i * Nb + i] * rmatCur[i]+rMat[i * Nb + i+1] * rmatCur[i+1]  + constArr[i]) + curRnd[i];
        }
        f1[Nb-1] = (ks / zeta) * (rMat[(Nb-1)*Nb+Nb-2] * rmatCur[Nb-2]+rMat[(Nb-1)*Nb+Nb-1] * rmatCur[Nb-1]  + constArr[Nb-1]) + curRnd[Nb-1];

//         printf("Initial rmatCur:\n");
//         for (int i = 0; i < Nb; i++) {
//             printf("%.8f ", f1[i]);
//         }
        // Predictor step (initial guess)
        for (int i = 0; i < Nb; i++) {
            initialGuess[i] = rmatCur[i] + h * f1[i];
        }

        // Corrector step (f2)
        f2[0] = (ks / zeta) * (rMat[0] * initialGuess[0] +rMat[1] * initialGuess[1]+ constArr[0]) + curRnd[0];

        for (int i = 1; i < Nb-1; i++) {
            f2[i] = (ks / zeta) * (rMat[i * Nb + i-1] * initialGuess[i-1]+rMat[i * Nb + i] * initialGuess[i] +rMat[i * Nb + i+1] * initialGuess[i+1]+ constArr[i]) + curRnd[i];
        }
        f2[Nb-1] = (ks / zeta) * (rMat[(Nb-1) * Nb + Nb-2] * initialGuess[Nb-2]+rMat[(Nb-1) * Nb + Nb-1] * initialGuess[Nb-1] + constArr[Nb-1]) + curRnd[Nb-1];

        // Update rmatNext using predictor/corrector
        for (int i = 0; i < Nb; i++) {
            rmatNext[i] = rmatCur[i] + 0.5 * h * (f1[i] + f2[i]);
        }

// 
//         // Frame selection: Store bead positions in sol array at selected frames
        if (k < numFrames-1 && t >= solFrames[k]) {
            // Store updated positions at subsequent frames
            for (int i = 0; i < Nb; i++) {
                sol[(k+1) * Nb + i] = rmatNext[i];
            }
            k++;
        }
// 
//         // Update rmatCur for the next iteration
        for (int i = 0; i < Nb; i++) {
            rmatCur[i] = rmatNext[i];
        }
//         memcpy(rmatCur, rmatNext, Nb * sizeof(double));
// 
//         // Add checks to prevent divergence (NaN or Inf)
        for (int i = 0; i < Nb; i++) {
            if (isnan(rmatCur[i]) || isinf(rmatCur[i])) {
                mexErrMsgIdAndTxt("MATLAB:rouse_sim:divergence", "Detected NaN or Inf value in rmatCur. Simulation diverged.");
            }
        }
    }

    // Free allocated memory
    free(rmatCur);
    free(rmatNext);
    free(rMat);
    free(constArr);
    free(curRnd);
    free(f1);
    free(f2);
    free(initialGuess);
    free(solFrames);
}

// Gateway function to interface with MATLAB
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Check for correct number of arguments
    if (nrhs != 9) {
        mexErrMsgIdAndTxt("MATLAB:rouse_sim:nrhs", "Nine input arguments required.");
    }

    // Extract input arguments
    double Xs = mxGetScalar(prhs[0]);
    int Nb = (int)mxGetScalar(prhs[1]);
    double zeta = mxGetScalar(prhs[2]);
    double ks = mxGetScalar(prhs[3]);
    double kT = mxGetScalar(prhs[4]);
    double h = mxGetScalar(prhs[5]);
    double fTime = mxGetScalar(prhs[6]);
    int frameLen = (int)mxGetScalar(prhs[7]);
    int numFrames = (int)mxGetScalar(prhs[8]);

    // Create an output matrix
    plhs[0] = mxCreateDoubleMatrix(Nb,numFrames, mxREAL);

    // Call the main simulation function
    double *sol = mxGetPr(plhs[0]);
    rouse_sim(Xs, Nb, zeta, ks, kT, h, fTime, frameLen, numFrames, sol);
}
