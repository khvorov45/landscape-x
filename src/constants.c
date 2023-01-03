#include "mltaln.h"
#include "miyata.h"
#include "miyata5.h"
#include "DNA.h"

#include "JTT.c"

#define DEFAULTGOP_B -1530
#define DEFAULTGEP_B -00
#define DEFAULTOFS_B -123

void
BLOSUMmtx(int n, double** matrix, double* freq, unsigned char* amino, char* amino_grp, int* rescalept) {
    char   locaminod[] = "ARNDCQEGHILKMFPSTWYVBZX.-J";
    char   locgrpd[] = {0, 3, 2, 2, 5, 2, 2, 0, 3, 1, 1, 3, 1, 4, 0, 0, 0, 4, 4, 1, 2, 2, 6, 6, 6, 1};
    double freqd[20] = {0.077, 0.051, 0.043, 0.052, 0.020, 0.041, 0.062, 0.074, 0.023, 0.052, 0.091, 0.059, 0.024, 0.040, 0.051, 0.069, 0.059, 0.014, 0.032, 0.066};

    // clang-format off
	double tmpmtx30[] = {
    4,
   -1,     8,
    0,    -2,     8,
    0,    -1,     1,     9,
   -3,    -2,    -1,    -3,    17,
    1,     3,    -1,    -1,    -2,     8,
    0,    -1,    -1,     1,     1,     2,     6,
    0,    -2,     0,    -1,    -4,    -2,    -2,     8,
   -2,    -1,    -1,    -2,    -5,     0,     0,    -3,    14,
    0,    -3,     0,    -4,    -2,    -2,    -3,    -1,    -2,     6,
   -1,    -2,    -2,    -1,     0,    -2,    -1,    -2,    -1,     2,     4,
    0,     1,     0,     0,    -3,     0,     2,    -1,    -2,    -2,    -2,     4,
    1,     0,     0,    -3,    -2,    -1,    -1,    -2,     2,     1,     2,     2,     6,
   -2,    -1,    -1,    -5,    -3,    -3,    -4,    -3,    -3,     0,     2,    -1,    -2,    10,
   -1,    -1,    -3,    -1,    -3,     0,     1,    -1,     1,    -3,    -3,     1,    -4,    -4,    11,
    1,    -1,     0,     0,    -2,    -1,     0,     0,    -1,    -1,    -2,     0,    -2,    -1,    -1,     4,
    1,    -3,     1,    -1,    -2,     0,    -2,    -2,    -2,     0,     0,    -1,     0,    -2,     0,     2,     5,
   -5,     0,    -7,    -4,    -2,    -1,    -1,     1,    -5,    -3,    -2,    -2,    -3,     1,    -3,    -3,    -5,    20,
   -4,     0,    -4,    -1,    -6,    -1,    -2,    -3,     0,    -1,     3,    -1,    -1,     3,    -2,    -2,    -1,     5,     9,
    1,    -1,    -2,    -2,    -2,    -3,    -3,    -3,    -3,     4,     1,    -2,     0,     1,    -4,    -1,     1,    -3,     1,     5,
    0,    -2,     4,     5,    -2,    -1,     0,     0,    -2,    -2,    -1,     0,    -2,    -3,    -2,     0,     0,    -5,    -3,    -2,     5,
    0,     0,    -1,     0,     0,     4,     5,    -2,     0,    -3,    -1,     1,    -1,    -4,     0,    -1,    -1,    -1,    -2,    -3,     0,     4,
    0,    -1,     0,    -1,    -2,     0,    -1,    -1,    -1,     0,     0,     0,     0,    -1,    -1,     0,     0,    -2,    -1,     0,    -1,     0,    -1};
	
	double tmpmtx45[] = {
    5,
   -2,      7,
   -1,      0,      6,
   -2,     -1,      2,      7,
   -1,     -3,     -2,     -3,     12,
   -1,      1,      0,      0,     -3,      6,
   -1,      0,      0,      2,     -3,      2,      6,
    0,     -2,      0,     -1,     -3,     -2,     -2,      7,
   -2,      0,      1,      0,     -3,      1,      0,     -2,     10,
   -1,     -3,     -2,     -4,     -3,     -2,     -3,     -4,     -3,      5,
   -1,     -2,     -3,     -3,     -2,     -2,     -2,     -3,     -2,      2,      5,
   -1,      3,      0,      0,     -3,      1,      1,     -2,     -1,     -3,     -3,      5,
   -1,     -1,     -2,     -3,     -2,      0,     -2,     -2,      0,      2,      2,     -1,      6,
   -2,     -2,     -2,     -4,     -2,     -4,     -3,     -3,     -2,      0,      1,     -3,      0,      8,
   -1,     -2,     -2,     -1,     -4,     -1,      0,     -2,     -2,     -2,     -3,     -1,     -2,     -3,      9,
    1,     -1,      1,      0,     -1,      0,      0,      0,     -1,     -2,     -3,     -1,     -2,     -2,     -1,      4,
    0,     -1,      0,     -1,     -1,     -1,     -1,     -2,     -2,     -1,     -1,     -1,     -1,     -1,     -1,      2,      5,
   -2,     -2,     -4,     -4,     -5,     -2,     -3,     -2,     -3,     -2,     -2,     -2,     -2,      1,     -3,     -4,     -3,     15,
   -2,     -1,     -2,     -2,     -3,     -1,     -2,     -3,      2,      0,      0,     -1,      0,      3,     -3,     -2,     -1,      3,      8,
    0,     -2,     -3,     -3,     -1,     -3,     -3,     -3,     -3,      3,      1,     -2,      1,      0,     -3,     -1,      0,     -3,     -1,      5};

    double tmpmtx50[] = {
    5,
   -2,      7,
   -1,     -1,      7,
   -2,     -2,      2,      8,
   -1,     -4,     -2,     -4,     13,
   -1,      1,      0,      0,     -3,      7,
   -1,      0,      0,      2,     -3,      2,      6,
    0,     -3,      0,     -1,     -3,     -2,     -3,      8,
   -2,      0,      1,     -1,     -3,      1,      0,     -2,     10,
   -1,     -4,     -3,     -4,     -2,     -3,     -4,     -4,     -4,      5,
   -2,     -3,     -4,     -4,     -2,     -2,     -3,     -4,     -3,      2,      5,
   -1,      3,      0,     -1,     -3,      2,      1,     -2,      0,     -3,     -3,      6,
   -1,     -2,     -2,     -4,     -2,      0,     -2,     -3,     -1,      2,      3,     -2,      7,
   -3,     -3,     -4,     -5,     -2,     -4,     -3,     -4,     -1,      0,      1,     -4,      0,      8,
   -1,     -3,     -2,     -1,     -4,     -1,     -1,     -2,     -2,     -3,     -4,     -1,     -3,     -4,     10,
    1,     -1,      1,      0,     -1,      0,     -1,      0,     -1,     -3,     -3,      0,     -2,     -3,     -1,      5,
    0,     -1,      0,     -1,     -1,     -1,     -1,     -2,     -2,     -1,     -1,     -1,     -1,     -2,     -1,      2,      5,
   -3,     -3,     -4,     -5,     -5,     -1,     -3,     -3,     -3,     -3,     -2,     -3,     -1,      1,     -4,     -4,     -3,     15,
   -2,     -1,     -2,     -3,     -3,     -1,     -2,     -3,      2,     -1,     -1,     -2,      0,      4,     -3,     -2,     -2,      2,      8,
    0,     -3,     -3,     -4,     -1,     -3,     -3,     -4,     -4,      4,      1,     -3,      1,     -1,     -3,     -2,      0,     -3,     -1,      5};

	double tmpmtx62[] = { 
    5.893685,
   -2.120252,  8.210189,
   -2.296072, -0.659672,  8.479856,
   -2.630151, -2.408668,  1.907550,  8.661363,
   -0.612761, -5.083814, -3.989626, -5.189966, 12.873172,
   -1.206025,  1.474162,  0.002529, -0.470069, -4.352838,  7.927704,
   -1.295821, -0.173087, -0.402015,  2.265459, -5.418729,  2.781955,  7.354247,
    0.239392, -3.456163, -0.634136, -1.970281, -3.750621, -2.677743, -3.165266,  8.344902,
   -2.437724, -0.374792,  0.867735, -1.678363, -4.481724,  0.672051, -0.176497, -3.061315, 11.266586,
   -1.982718, -4.485360, -4.825558, -4.681732, -1.841495, -4.154454, -4.791538, -5.587336, -4.847345,  5.997760,
   -2.196882, -3.231860, -5.068375, -5.408471, -1.916207, -3.200863, -4.269723, -5.440437, -4.180099,  2.282412,  5.774148,
   -1.101017,  3.163105, -0.268534, -1.052724, -4.554510,  1.908859,  1.163010, -2.291924, -1.081539, -4.005209, -3.670219,  6.756827,
   -1.402897, -2.050705, -3.226290, -4.587785, -2.129758, -0.631437, -2.997038, -4.014898, -2.326896,  1.690191,  2.987638, -2.032119,  8.088951,
   -3.315080, -4.179521, -4.491005, -5.225795, -3.563219, -4.746598, -4.788639, -4.661029, -1.851231, -0.241317,  0.622170, -4.618016,  0.018880,  9.069126,
   -1.221394, -3.162863, -3.000581, -2.220163, -4.192770, -1.922917, -1.674258, -3.200320, -3.241363, -4.135001, -4.290107, -1.520445, -3.714633, -5.395930, 11.046892,
    1.673639, -1.147170,  0.901353, -0.391548, -1.312485, -0.151708, -0.220375, -0.438748, -1.322366, -3.522266, -3.663923, -0.305170, -2.221304, -3.553533, -1.213470,  5.826527,
   -0.068042, -1.683495, -0.069138, -1.576054, -1.299983, -1.012997, -1.294878, -2.363065, -2.528844, -1.076382, -1.796229, -1.004336, -0.999449, -3.161436, -1.612919,  2.071710,  6.817956,
   -3.790328, -4.019108, -5.543911, -6.321502, -3.456164, -2.919725, -4.253197, -3.737232, -3.513238, -3.870811, -2.447829, -4.434676, -2.137255,  1.376341, -5.481260, -4.127804, -3.643382, 15.756041,
   -2.646022, -2.540799, -3.122641, -4.597428, -3.610671, -2.131601, -3.030688, -4.559647,  2.538948, -1.997058, -1.593097, -2.730047, -1.492308,  4.408690, -4.379667, -2.528713, -2.408996,  3.231335,  9.892544,
   -0.284140, -3.753871, -4.314525, -4.713963, -1.211518, -3.297575, -3.663425, -4.708118, -4.676220,  3.820569,  1.182672, -3.393535,  1.030861, -1.273542, -3.523054, -2.469318, -0.083276, -4.251392, -1.811267,  5.653391};

	double tmpmtx80[] = {
    7,
   -3,      9,
   -3,     -1,      9,
   -3,     -3,      2,     10,
   -1,     -6,     -5,     -7,     13,
   -2,      1,      0,     -1,     -5,      9,
   -2,     -1,     -1,      2,     -7,      3,      8,
    0,     -4,     -1,     -3,     -6,     -4,     -4,      9,
   -3,      0,      1,     -2,     -7,      1,      0,     -4,     12,
   -3,     -5,     -6,     -7,     -2,     -5,     -6,     -7,     -6,      7,
   -3,     -4,     -6,     -7,     -3,     -4,     -6,     -7,     -5,      2,      6,
   -1,      3,      0,     -2,     -6,      2,      1,     -3,     -1,     -5,     -4,      8,
   -2,     -3,     -4,     -6,     -3,     -1,     -4,     -5,     -4,      2,      3,     -3,      9,
   -4,     -5,     -6,     -6,     -4,     -5,     -6,     -6,     -2,     -1,      0,     -5,      0,     10,
   -1,     -3,     -4,     -3,     -6,     -3,     -2,     -5,     -4,     -5,     -5,     -2,     -4,     -6,     12,
    2,     -2,      1,     -1,     -2,     -1,     -1,     -1,     -2,     -4,     -4,     -1,     -3,     -4,     -2,      7,
    0,     -2,      0,     -2,     -2,     -1,     -2,     -3,     -3,     -2,     -3,     -1,     -1,     -4,     -3,      2,      8,
   -5,     -5,     -7,     -8,     -5,     -4,     -6,     -6,     -4,     -5,     -4,     -6,     -3,      0,     -7,     -6,     -5,     16,
   -4,     -4,     -4,     -6,     -5,     -3,     -5,     -6,      3,     -3,     -2,     -4,     -3,      4,     -6,     -3,     -3,      3,     11,
   -1,     -4,     -5,     -6,     -2,     -4,     -4,     -6,     -5,      4,      1,     -4,      1,     -2,     -4,     -3,      0,     -5,     -3,      7};

	double tmpmtx90[] = {
    5,
   -2,  6,
   -2, -1,  7,
   -3, -3,  1,  7,
   -1, -5, -4, -5,  9,
   -1,  1,  0, -1, -4,  7,
   -1, -1, -1,  1, -6,  2,  6,
    0, -3, -1, -2, -4, -3, -3,  6,
   -2,  0,  0, -2, -5,  1, -1, -3,  8,
   -2, -4, -4, -5, -2, -4, -4, -5, -4,  5,
   -2, -3, -4, -5, -2, -3, -4, -5, -4,  1,  5,
   -1,  2,  0, -1, -4,  1,  0, -2, -1, -4, -3,  6,
   -2, -2, -3, -4, -2,  0, -3, -4, -3,  1,  2, -2,  7,
   -3, -4, -4, -5, -3, -4, -5, -5, -2, -1,  0, -4, -1,  7,
   -1, -3, -3, -3, -4, -2, -2, -3, -3, -4, -4, -2, -3, -4,  8,
    1, -1,  0, -1, -2, -1, -1, -1, -2, -3, -3, -1, -2, -3, -2,  5,
    0, -2,  0, -2, -2, -1, -1, -3, -2, -1, -2, -1, -1, -3, -2,  1,  6,
   -4, -4, -5, -6, -4, -3, -5, -4, -3, -4, -3, -5, -2,  0, -5, -4, -4, 11,
   -3, -3, -3, -4, -4, -3, -4, -5,  1, -2, -2, -3, -2,  3, -4, -3, -2,  2,  8,
   -1, -3, -4, -5, -2, -3, -3, -5, -4,  3,  0, -3,  0, -2, -3, -2, -1, -3, -3,  5};

	double tmpmtx100[] = {
    8,
   -3,10,
   -4,-2,11,
   -5,-5, 1,10,
   -2,-8,-5,-8,14,
   -2, 0,-1,-2,-7,11,
   -3,-2,-2, 2,-9, 2,10,
   -1,-6,-2,-4,-7,-5,-6, 9,
   -4,-1, 0,-3,-8, 1,-2,-6,13,
   -4,-7,-7,-8,-3,-6,-7,-9,-7, 8,
   -4,-6,-7,-8,-5,-5,-7,-8,-6, 2, 8,
   -2, 3,-1,-3,-8, 2, 0,-5,-3,-6,-6,10,
   -3,-4,-5,-8,-4,-2,-5,-7,-5, 1, 3,-4,12,
   -5,-6,-7,-8,-4,-6,-8,-8,-4,-2, 0,-6,-1,11,
   -2,-5,-5,-5,-8,-4,-4,-6,-5,-7,-7,-3,-5,-7,12,
    1,-3, 0,-2,-3,-2,-2,-2,-3,-5,-6,-2,-4,-5,-3, 9,
   -1,-3,-1,-4,-3,-3,-3,-5,-4,-3,-4,-3,-2,-5,-4, 2, 9,
   -6,-7,-8,-10,-7,-5,-8,-7,-5,-6,-5,-8,-4, 0,-8,-7,-7,17,
   -5,-5,-5,-7,-6,-4,-7,-8, 1,-4,-4,-5,-5, 4,-7,-5,-5, 2,12,
   -2,-6,-7,-8,-3,-5,-5,-8,-7, 4, 0,-5, 0,-3,-6,-4,-1,-5,-5, 8};

	double tmpmtx0[] = {
    2.4,
   -0.6,    4.7,
   -0.3,    0.3,    3.8,
   -0.3,   -0.3,    2.2,    4.7,
    0.5,   -2.2,   -1.8,   -3.2,   11.5,
   -0.2,    1.5,    0.7,    0.9,   -2.4,    2.7,
    0.0,    0.4,    0.9,    2.7,   -3.0,    1.7,    3.6,
    0.5,   -1.0,    0.4,    0.1,   -2.0,   -1.0,   -0.8,    6.6,
   -0.8,    0.6,    1.2,    0.4,   -1.3,    1.2,    0.4,   -1.4,    6.0,
   -0.8,   -2.4,   -2.8,   -3.8,   -1.1,   -1.9,   -2.7,   -4.5,   -2.2,    4.0,
   -1.2,   -2.2,   -3.0,   -4.0,   -1.5,   -1.6,   -2.8,   -4.4,   -1.9,    2.8,    4.0,
   -0.4,    2.7,    0.8,    0.5,   -2.8,    1.5,    1.2,   -1.1,    0.6,   -2.1,   -2.1,    3.2,
   -0.7,   -1.7,   -2.2,   -3.0,   -0.9,   -1.0,   -2.0,   -3.5,   -1.3,    2.5,    2.8,   -1.4,    4.3,
   -2.3,   -3.2,   -3.1,   -4.5,   -0.8,   -2.6,   -3.9,   -5.2,   -0.1,    1.0,    2.0,   -3.3,    1.6,    7.0,
    0.3,   -0.9,   -0.9,   -0.7,   -3.1,   -0.2,   -0.5,   -1.6,   -1.1,   -2.6,   -2.3,   -0.6,   -2.4,   -3.8,    7.6,
    1.1,   -0.2,    0.9,    0.5,    0.1,    0.2,    0.2,    0.4,   -0.2,   -1.8,   -2.1,    0.1,   -1.4,   -2.8,    0.4,    2.2,
    0.6,   -0.2,    0.5,    0.0,   -0.5,    0.0,   -0.1,   -1.1,   -0.3,   -0.6,   -1.3,    0.1,   -0.6,   -2.2,    0.1,    1.5,    2.5,
   -3.6,   -1.6,   -3.6,   -5.2,   -1.0,   -2.7,   -4.3,   -4.0,   -0.8,   -1.8,   -0.7,   -3.5,   -1.0,    3.6,   -5.0,   -3.3,   -3.5,   14.2,
   -2.2,   -1.8,   -1.4,   -2.8,   -0.5,   -1.7,   -2.7,   -4.0,    2.2,   -0.7,    0.0,   -2.1,   -0.2,    5.1,   -3.1,   -1.9,   -1.9,    4.1,    7.8,
    0.1,   -2.0,   -2.2,   -2.9,    0.0,   -1.5,   -1.9,   -3.3,   -2.0,    3.1,    1.8,   -1.7,    1.6,    0.1,   -1.8,   -1.0,    0.0,   -2.6,   -1.1,    3.4};
    // clang-format on

    int     i, j, count;
    double  av;
    double* tmpmtx;

    if (n == 30)
        tmpmtx = tmpmtx30;
    else if (n == 45)
        tmpmtx = tmpmtx45;
    else if (n == 50)
        tmpmtx = tmpmtx50;
    else if (n == 62)
        tmpmtx = tmpmtx62;
    else if (n == 80)
        tmpmtx = tmpmtx80;
    else if (n == 90)
        tmpmtx = tmpmtx90;
    else if (n == 100)
        tmpmtx = tmpmtx100;
    else if (n == 0)
        tmpmtx = tmpmtx0;
    else if (n == -1)
        tmpmtx = loadaamtx(rescalept);
    else {
        fprintf(stderr, "blosum %d ?\n", n);
        exit(1);
    }

    count = 0;
    for (i = 0; i < 20; i++) {
        for (j = 0; j <= i; j++) {
            matrix[i][j] = matrix[j][i] = (double)tmpmtx[count++];
        }
    }
    if (n == -1 && tmpmtx[400] != -1.0) {
        for (i = 0; i < 20; i++)
            freq[i] = tmpmtx[400 + i];
        av = 0.0;
        for (i = 0; i < 20; i++)
            av += freq[i];
        for (i = 0; i < 20; i++)
            freq[i] /= av;
    } else
        for (i = 0; i < 20; i++)
            freq[i] = freqd[i];

    if (n == -1)
        free(tmpmtx);

    for (i = 0; i < 26; i++)
        amino[i] = locaminod[i];
    for (i = 0; i < 26; i++)
        amino_grp[(int)amino[i]] = locgrpd[i];
}

static int
checkchar(int i) {
    if (i > 0xff || i <= 0x00 || i == 0x3E || i == 0x3D || i == 0x3C || i == 0x2D || i == 0x20 || i == 0x0d || i == 0x0a)
        return (1);
    return (0);
}

static int
overridematrix(double** matrix) {
    char         buf[500];
    FILE*        fp;
    unsigned int i1, i2;
    double       v;
    char*        bpt;
    int          nread;
    int          userdefined;

    fp = fopen("_aamtx", "r");
    if (fp == NULL) {
        fprintf(stderr, "warning: cannot open scorematrix. Use the default one.\n");
        //		f2cl.c de tomaranai youni
        //		exit( 1 );
        return 0;
    }

    userdefined = 0;
    while (1) {
        fgets(buf, 499, fp);
        if (feof(fp))
            break;

        if ((bpt = strchr(buf, '#')))
            *bpt = 0;

        i1 = i2 = 0;
        nread = sscanf(buf, "%x %x %lf", &i1, &i2, &v);
        if (nread == EOF)
            continue;
        if (nread != 3) {
            reporterr("Format error in this line?\n");
            reporterr("%s\n", buf);
            //			reporterr( "To set a score, 100, for a match of 0x41 and 0x42, \n" );
            //			reporterr( "0x41 0x42 100 \n" );
            exit(1);
        }
        if (checkchar(i1)) {
            reporterr("%c=0x%x cannot be used (1)\n", i1, i1);
            exit(1);
        }
        if (checkchar(i2)) {
            reporterr("%c=0x%x cannot be used (2)\n", i2, i2);
            exit(1);
        }
        //		reporterr( "Score(%c=0x%x,%c=0x%x)=%f\n", i1, i1, i2, i2, v );
        matrix[i1][i2] = v;
        matrix[i2][i1] = v;  // 2018/May/11
        userdefined = 1;
    }
    fclose(fp);
    return userdefined;
}

int
extendedmtx(double** matrix, double* freq, unsigned char* amino, char* amino_grp) {
    int i;
    int j;
    int userdefined;

    for (i = 0; i < nalphabets; i++) {
        //		fprintf( stderr, "i=%d, i=%c\n", i, i );
        amino[i] = (unsigned char)i;
    }
    for (i = 0; i < nalphabets; i++)
        amino_grp[(int)amino[i]] = i % 6;
    for (i = 0; i < nalphabets; i++)
        freq[i] = 1.0 / nalphabets;

    for (i = 0; i < nalphabets; i++) {
        for (j = 0; j <= i; j++) {
            matrix[i][j] = matrix[j][i] = (double)-1.0;
        }
    }
    //	for( i=0; i<nalphabets; i++ )
    //		matrix[i][i] = matrix[i][i] = (double)1.0;

    userdefined = overridematrix(matrix);

#if 0  // user-defined matrix + user-defined freqency ga aru toki fukkatsu saseru.
	if( tmpmtx[400] != -1.0 ) 
	{
		for( i=0; i<20; i++ ) freq[i] = tmpmtx[400+i];
		av = 0.0;
		for( i=0; i<20; i++ ) av += freq[i];
		for( i=0; i<20; i++ ) freq[i] /= av;
	}
	else
		for( i=0; i<20; i++ ) freq[i] = freqd[i];
#endif
#if 0
	for( i=0; i<nalphabets; i++ )
	{
		fprintf( stderr, "%d: %c, %d, %f\n", i, amino[i], amino_grp[amino[i]], freq[i] );
	}
#endif
    return userdefined;
}

#define DEBUG 0
#define TEST 0

#define NORMALIZE1 1

static int
shishagonyuu(double in) {
    int out;
    if (in > 0.0)
        out = ((int)(in + 0.5));
    else if (in == 0.0)
        out = (0);
    else if (in < 0.0)
        out = ((int)(in - 0.5));
    else
        out = 0;
    return (out);
}

static void
nscore(int* amino_n, int** n_dis) {
    int i;
    for (i = 0; i < 26; i++) {
        //		reporterr( "i=%d (%c), n_dis[%d][%d] = %d\n", i, amino[i], i, amino_n['n'], n_dis[i][amino_n['n']] );
        n_dis[i][amino_n['n']] = shishagonyuu((double)0.25 * n_dis[i][i]);
        //		reporterr( "-> i=%d, n_dis[%d][%d] = %d\n", i, i, amino_n['n'], n_dis[i][amino_n['n']] );
        n_dis[amino_n['n']][i] = n_dis[i][amino_n['n']];
    }
    //	n_dis[amino_n['n']][amino_n['n']] = shishagonyuu( (double)0.25 * 0.25 * ( n_dis[0][0] + n_dis[1][1] + n_dis[2][2] + n_dis[3][3] ) );
    n_dis[amino_n['n']][amino_n['n']] = shishagonyuu((double)0.25 * (n_dis[0][0] + n_dis[1][1] + n_dis[2][2] + n_dis[3][3]));  // 2017/Jan/2

#if 0  // Ato de kakunin
	for( i=0; i<26; i++ )
	{
		n_dis[i][amino_n['-']] = shishagonyuu( (double)0.25 * n_dis[i][i] );
		n_dis[amino_n['-']][i] = n_dis[i][amino_n['-']];
	}
//	n_dis[amino_n['-']][amino_n['-']] = shishagonyuu( (double)0.25 * 0.25 * ( n_dis[0][0] + n_dis[1][1] + n_dis[2][2] + n_dis[3][3] ) ); // DAME!
#endif
}

static void
ambiguousscore(int* amino_n, int** n_dis) {
    int i;
    for (i = 0; i < 26; i++) {
        n_dis[i][amino_n['r']] = shishagonyuu((double)1 / 2 * (n_dis[amino_n['a']][i] + n_dis[amino_n['g']][i]));
        n_dis[i][amino_n['y']] = shishagonyuu((double)1 / 2 * (n_dis[amino_n['c']][i] + n_dis[amino_n['t']][i]));
        n_dis[i][amino_n['k']] = shishagonyuu((double)1 / 2 * (n_dis[amino_n['g']][i] + n_dis[amino_n['t']][i]));
        n_dis[i][amino_n['m']] = shishagonyuu((double)1 / 2 * (n_dis[amino_n['a']][i] + n_dis[amino_n['c']][i]));
        n_dis[i][amino_n['s']] = shishagonyuu((double)1 / 2 * (n_dis[amino_n['g']][i] + n_dis[amino_n['c']][i]));
        n_dis[i][amino_n['w']] = shishagonyuu((double)1 / 2 * (n_dis[amino_n['a']][i] + n_dis[amino_n['t']][i]));
        n_dis[i][amino_n['b']] = shishagonyuu((double)1 / 3 * (n_dis[amino_n['c']][i] + n_dis[amino_n['g']][i] + n_dis[amino_n['t']][i]));
        n_dis[i][amino_n['d']] = shishagonyuu((double)1 / 3 * (n_dis[amino_n['a']][i] + n_dis[amino_n['g']][i] + n_dis[amino_n['t']][i]));
        n_dis[i][amino_n['h']] = shishagonyuu((double)1 / 3 * (n_dis[amino_n['a']][i] + n_dis[amino_n['c']][i] + n_dis[amino_n['t']][i]));
        n_dis[i][amino_n['v']] = shishagonyuu((double)1 / 3 * (n_dis[amino_n['a']][i] + n_dis[amino_n['c']][i] + n_dis[amino_n['g']][i]));

        n_dis[amino_n['r']][i] = n_dis[i][amino_n['r']];
        n_dis[amino_n['y']][i] = n_dis[i][amino_n['y']];
        n_dis[amino_n['k']][i] = n_dis[i][amino_n['k']];
        n_dis[amino_n['m']][i] = n_dis[i][amino_n['m']];
        n_dis[amino_n['s']][i] = n_dis[i][amino_n['s']];
        n_dis[amino_n['w']][i] = n_dis[i][amino_n['w']];
        n_dis[amino_n['b']][i] = n_dis[i][amino_n['b']];
        n_dis[amino_n['d']][i] = n_dis[i][amino_n['d']];
        n_dis[amino_n['h']][i] = n_dis[i][amino_n['h']];
        n_dis[amino_n['v']][i] = n_dis[i][amino_n['v']];
    }

    i = amino_n['r'];
    n_dis[i][i] = shishagonyuu((double)1 / 2 * (n_dis[amino_n['a']][amino_n['a']] + n_dis[amino_n['g']][amino_n['g']]));
    i = amino_n['y'];
    n_dis[i][i] = shishagonyuu((double)1 / 2 * (n_dis[amino_n['c']][amino_n['c']] + n_dis[amino_n['t']][amino_n['t']]));
    i = amino_n['k'];
    n_dis[i][i] = shishagonyuu((double)1 / 2 * (n_dis[amino_n['g']][amino_n['g']] + n_dis[amino_n['t']][amino_n['t']]));
    i = amino_n['m'];
    n_dis[i][i] = shishagonyuu((double)1 / 2 * (n_dis[amino_n['a']][amino_n['a']] + n_dis[amino_n['c']][amino_n['c']]));
    i = amino_n['s'];
    n_dis[i][i] = shishagonyuu((double)1 / 2 * (n_dis[amino_n['g']][amino_n['g']] + n_dis[amino_n['c']][amino_n['c']]));
    i = amino_n['w'];
    n_dis[i][i] = shishagonyuu((double)1 / 2 * (n_dis[amino_n['a']][amino_n['a']] + n_dis[amino_n['t']][amino_n['t']]));
    i = amino_n['b'];
    n_dis[i][i] = shishagonyuu((double)1 / 3 * (n_dis[amino_n['c']][amino_n['c']] + n_dis[amino_n['g']][amino_n['g']] + n_dis[amino_n['t']][amino_n['t']]));
    i = amino_n['d'];
    n_dis[i][i] = shishagonyuu((double)1 / 3 * (n_dis[amino_n['a']][amino_n['a']] + n_dis[amino_n['g']][amino_n['g']] + n_dis[amino_n['t']][amino_n['t']]));
    i = amino_n['h'];
    n_dis[i][i] = shishagonyuu((double)1 / 3 * (n_dis[amino_n['a']][amino_n['a']] + n_dis[amino_n['c']][amino_n['c']] + n_dis[amino_n['t']][amino_n['t']]));
    i = amino_n['v'];
    n_dis[i][i] = shishagonyuu((double)1 / 3 * (n_dis[amino_n['a']][amino_n['a']] + n_dis[amino_n['c']][amino_n['c']] + n_dis[amino_n['g']][amino_n['g']]));
}

static void
calcfreq_nuc(int nseq, char** seq, double* datafreq) {
    int    i, j, l;
    int    aan;
    double total;
    for (i = 0; i < 4; i++)
        datafreq[i] = 0.0;
    total = 0.0;
    for (i = 0; i < nseq; i++) {
        l = strlen(seq[i]);
        for (j = 0; j < l; j++) {
            aan = amino_n[(int)seq[i][j]];
            if (aan == 4)
                aan = 3;
            if (aan >= 0 && aan < 4) {
                datafreq[aan] += 1.0;
                total += 1.0;
            }
        }
    }
    total = 0.0;
    for (i = 0; i < 4; i++)
        total += datafreq[i];
    for (i = 0; i < 4; i++)
        datafreq[i] /= (double)total;
    for (i = 0; i < 4; i++)
        if (datafreq[i] < 0.0001)
            datafreq[i] = 0.0001;

    total = 0.0;
    for (i = 0; i < 4; i++)
        total += datafreq[i];
    //	reporterr(       "total = %f\n", total );
    for (i = 0; i < 4; i++)
        datafreq[i] /= (double)total;

#if 0
	reporterr(       "\ndatafreq = " );
	for( i=0; i<4; i++ )
		reporterr(       "%10.5f ", datafreq[i] );
	reporterr(       "\n" );
	exit( 1 );
#endif
}

static void
calcfreq(int nseq, char** seq, double* datafreq) {
    int    i, j, l;
    int    aan;
    double total;
    for (i = 0; i < nscoredalphabets; i++)
        datafreq[i] = 0.0;
    total = 0.0;
    for (i = 0; i < nseq; i++) {
        l = strlen(seq[i]);
        for (j = 0; j < l; j++) {
            aan = amino_n[(int)seq[i][j]];
            if (aan >= 0 && aan < nscoredalphabets && seq[i][j] != '-') {
                datafreq[aan] += 1.0;
                total += 1.0;
            }
        }
    }
    total = 0.0;
    for (i = 0; i < nscoredalphabets; i++)
        total += datafreq[i];
    for (i = 0; i < nscoredalphabets; i++)
        datafreq[i] /= (double)total;
    for (i = 0; i < nscoredalphabets; i++)
        if (datafreq[i] < 0.0001)
            datafreq[i] = 0.0001;

    //	reporterr(       "datafreq = \n" );
    //	for( i=0; i<nscoredalphabets; i++ )
    //		reporterr(       "%f\n", datafreq[i] );

    total = 0.0;
    for (i = 0; i < nscoredalphabets; i++)
        total += datafreq[i];
    //	reporterr(       "total = %f\n", total );
    for (i = 0; i < nscoredalphabets; i++)
        datafreq[i] /= (double)total;
}

void
calcfreq_from_scoremtx(double** n_distmp, double* datafreq) {
    int i, j;
    int nused = 0;
    for (i = 0; i < nalphabets; i++)
        datafreq[i] = 0.0;
    for (i = 0; i < nalphabets; i++)
        for (j = 0; j < i; j++) {
            if (n_distmp[i][j] != -1) {
                if (datafreq[i] == 0.0)
                    nused += 1;
                if (datafreq[j] == 0.0)
                    nused += 1;
                datafreq[i] = datafreq[j] = 1.0;
            }
        }
    for (i = 0; i < nalphabets; i++)
        datafreq[i] /= (double)nused;
    reporterr("nused=\n", nused);
    //	for( i=0; i<nalphabets; i++ ) reporterr( "%f\n", datafreq[i] );
}

static int
checkscoremtx(double** n_distmp, int nseq, char** seq) {
    int i, j, l, k;
    int aan;
    for (i = 0; i < nseq; i++) {
        l = strlen(seq[i]);
        for (j = 0; j < l; j++) {
            aan = amino_n[(unsigned char)seq[i][j]];
            for (k = 0; k < nalphabets; k++) {
                //				if( n_distmp[k][aan] != -1.0 && k != aan ) break;
                if (n_distmp[k][aan] != -1.0)
                    break;
            }
            if (k == nalphabets)
                return aan;
        }
    }
    return 0;
}

static void
calcfreq_extended(int nseq, char** seq, double* datafreq) {
    int    i, j, l;
    int    aan;
    double total;
    for (i = 0; i < nscoredalphabets; i++)
        datafreq[i] = 0.0;
    total = 0.0;
    for (i = 0; i < nseq; i++) {
        l = strlen(seq[i]);
        for (j = 0; j < l; j++) {
            aan = amino_n[(unsigned char)seq[i][j]];
            if (aan >= 0 && aan < nscoredalphabets && seq[i][j] != '-') {
                datafreq[aan] += 1.0;
                total += 1.0;
            }
        }
    }
    total = 0.0;
    for (i = 0; i < nscoredalphabets; i++)
        total += datafreq[i];
    for (i = 0; i < nscoredalphabets; i++)
        datafreq[i] /= (double)total;
        //	for( i=0; i<nscoredalphabets; i++ ) if( datafreq[i] < 0.0001 ) datafreq[i] = 0.0001;

#if 0
	reporterr(       "datafreq = \n" );
	for( i=0; i<nscoredalphabets; i++ )
		reporterr(       "%d %c %f\n", i, amino[i], datafreq[i] );
#endif

    total = 0.0;
    for (i = 0; i < nscoredalphabets; i++)
        total += datafreq[i];
    //	reporterr(       "total = %f\n", total );
    for (i = 0; i < nscoredalphabets; i++)
        datafreq[i] /= (double)total;
}

static void
generatenuc1pam(double** pam1, int kimuraR, double* freq) {
    int    i, j;
    double R[4][4], mut[4], total, tmp;

    R[0][0] = 0.0;
    R[0][1] = kimuraR;
    R[0][2] = 1.0;
    R[0][3] = 1.0;
    R[1][0] = kimuraR;
    R[1][1] = 0.0;
    R[1][2] = 1.0;
    R[1][3] = 1.0;
    R[2][0] = 1.0;
    R[2][1] = 1.0;
    R[2][2] = 0.0;
    R[2][3] = kimuraR;
    R[3][0] = 1.0;
    R[3][1] = 1.0;
    R[3][2] = kimuraR;
    R[3][3] = 0.0;

    total = 0.0;
    for (i = 0; i < 4; i++) {
        tmp = 0.0;
        for (j = 0; j < 4; j++)
            tmp += R[i][j] * freq[j];
        mut[i] = tmp;
        total += tmp * freq[i];
    }
    for (i = 0; i < 4; i++)
        for (j = 0; j < 4; j++) {
            if (i != j)
                pam1[i][j] = 0.01 / total * R[i][j] * freq[j];
            else
                pam1[i][j] = 1.0 - 0.01 / total * mut[i];
        }
}

void
constants(int nseq, char** seq) {
    int i, j, x;
    //	double tmp;
    char shiftmodel[100];
    int  charsize;

    if (nblosum < 0)
        dorp = 'p';

    if (penalty_shift_factor >= 10)
        trywarp = 0;
    else
        trywarp = 1;

    if (dorp == 'd') /* DNA */
    {
        int      k, m;
        double   average;
        double** pamx = AllocateDoubleMtx(11, 11);
        double** pam1 = AllocateDoubleMtx(4, 4);
        double*  freq = AllocateDoubleVec(4);

        nalphabets = 26;
        nscoredalphabets = 10;
        charsize = 0x80;

        n_dis = AllocateIntMtx(nalphabets, nalphabets);
        n_disLN = AllocateDoubleMtx(nalphabets, nalphabets);

        scoremtx = -1;
        if (RNAppenalty == NOTSPECIFIED)
            RNAppenalty = DEFAULTRNAGOP_N;
        if (RNAppenalty_ex == NOTSPECIFIED)
            RNAppenalty_ex = DEFAULTRNAGEP_N;
        if (ppenalty == NOTSPECIFIED)
            ppenalty = DEFAULTGOP_N;
        if (ppenalty_dist == NOTSPECIFIED)
            ppenalty_dist = ppenalty;
        if (ppenalty_OP == NOTSPECIFIED)
            ppenalty_OP = DEFAULTGOP_N;
        if (ppenalty_ex == NOTSPECIFIED)
            ppenalty_ex = DEFAULTGEP_N;
        if (ppenalty_EX == NOTSPECIFIED)
            ppenalty_EX = DEFAULTGEP_N;
        if (poffset == NOTSPECIFIED)
            poffset = DEFAULTOFS_N;
        if (RNApthr == NOTSPECIFIED)
            RNApthr = DEFAULTRNATHR_N;
        if (pamN == NOTSPECIFIED)
            pamN = DEFAULTPAMN;
        if (kimuraR == NOTSPECIFIED)
            kimuraR = 2;

        RNApenalty = (int)(3 * 600.0 / 1000.0 * RNAppenalty + 0.5);
        RNApenalty_ex = (int)(3 * 600.0 / 1000.0 * RNAppenalty_ex + 0.5);
        //		reporterr(       "DEFAULTRNAGOP_N = %d\n", DEFAULTRNAGOP_N );
        //		reporterr(       "RNAppenalty = %d\n", RNAppenalty );
        //		reporterr(       "RNApenalty = %d\n", RNApenalty );

        RNAthr = (int)(3 * 600.0 / 1000.0 * RNApthr + 0.5);
        penalty = (int)(3 * 600.0 / 1000.0 * ppenalty + 0.5);
        penalty_dist = (int)(3 * 600.0 / 1000.0 * ppenalty_dist + 0.5);
        penalty_shift = (int)(penalty_shift_factor * penalty);
        penalty_OP = (int)(3 * 600.0 / 1000.0 * ppenalty_OP + 0.5);
        penalty_ex = (int)(3 * 600.0 / 1000.0 * ppenalty_ex + 0.5);
        penalty_EX = (int)(3 * 600.0 / 1000.0 * ppenalty_EX + 0.5);
        offset = (int)(1 * 600.0 / 1000.0 * poffset + 0.5);
        offsetFFT = (int)(1 * 600.0 / 1000.0 * (-0) + 0.5);
        offsetLN = (int)(1 * 600.0 / 1000.0 * 100 + 0.5);
        penaltyLN = (int)(3 * 600.0 / 1000.0 * -2000 + 0.5);
        penalty_exLN = (int)(3 * 600.0 / 1000.0 * -100 + 0.5);

        if (trywarp)
            sprintf(shiftmodel, "%4.2f (%4.2f)", -(double)penalty_shift / 1800, -(double)penalty_shift / 600);
        else
            sprintf(shiftmodel, "noshift");

        sprintf(modelname, "%s%d (%d), %4.2f (%4.2f), %4.2f (%4.2f), %s", rnakozo ? "RNA" : "DNA", pamN, kimuraR, -(double)ppenalty * 0.001, -(double)ppenalty * 0.003, -(double)poffset * 0.001, -(double)poffset * 0.003, shiftmodel);

        for (i = 0; i < 26; i++)
            amino[i] = locaminon[i];
        for (i = 0; i < 0x80; i++)
            amino_n[i] = -1;
        for (i = 0; i < 26; i++)
            amino_n[(int)amino[i]] = i;
        if (fmodel == 1) {
            calcfreq_nuc(nseq, seq, freq);
            reporterr("a, freq[0] = %f\n", freq[0]);
            reporterr("g, freq[1] = %f\n", freq[1]);
            reporterr("c, freq[2] = %f\n", freq[2]);
            reporterr("t, freq[3] = %f\n", freq[3]);
        } else {
            freq[0] = 0.25;
            freq[1] = 0.25;
            freq[2] = 0.25;
            freq[3] = 0.25;
        }

        if (kimuraR == 9999) {
            for (i = 0; i < 4; i++)
                for (j = 0; j < 4; j++)
                    pamx[i][j] = (double)locn_disn[i][j];
#if NORMALIZE1
            average = 0.0;
            for (i = 0; i < 4; i++)
                for (j = 0; j < 4; j++)
                    average += pamx[i][j];
            average /= 16.0;

            if (disp)
                reporterr("average = %f\n", average);

            for (i = 0; i < 4; i++)
                for (j = 0; j < 4; j++)
                    pamx[i][j] -= average;

            for (i = 0; i < 4; i++)
                for (j = 0; j < 4; j++)
                    pamx[i][j] *= 600.0 / average;

            for (i = 0; i < 4; i++)
                for (j = 0; j < 4; j++)
                    pamx[i][j] -= offset;
#endif
        } else {
#if 0
				double f = 0.99;
				double s = (double)kimuraR / ( 2 + kimuraR ) * 0.01;
				double v = (double)1       / ( 2 + kimuraR ) * 0.01;
				pam1[0][0] = f; pam1[0][1] = s; pam1[0][2] = v; pam1[0][3] = v;
				pam1[1][0] = s; pam1[1][1] = f; pam1[1][2] = v; pam1[1][3] = v;
				pam1[2][0] = v; pam1[2][1] = v; pam1[2][2] = f; pam1[2][3] = s;
				pam1[3][0] = v; pam1[3][1] = v; pam1[3][2] = s; pam1[3][3] = f;
#else
            generatenuc1pam(pam1, kimuraR, freq);
#endif

            reporterr("generating a scoring matrix for nucleotide (dist=%d) ... ", pamN);

            if (disp) {
                reporterr(" TPM \n");
                for (i = 0; i < 4; i++) {
                    for (j = 0; j < 4; j++)
                        reporterr("%+#6.10f", pam1[i][j]);
                    reporterr("\n");
                }
                reporterr("\n");
            }

            MtxuntDouble(pamx, 4);
            for (x = 0; x < pamN; x++)
                MtxmltDouble(pamx, pam1, 4);

            if (disp) {
                reporterr(" TPM \n");
                for (i = 0; i < 4; i++) {
                    for (j = 0; j < 4; j++)
                        reporterr("%+#6.10f", pamx[i][j]);
                    reporterr("\n");
                }
                reporterr("\n");
            }

            for (i = 0; i < 4; i++)
                for (j = 0; j < 4; j++)
                    pamx[i][j] /= freq[j];
            //					pamx[i][j] /= 0.25;

            for (i = 0; i < 4; i++)
                for (j = 0; j < 4; j++) {
                    if (pamx[i][j] == 0.0) {
                        reporterr("WARNING: pamx[i][j] = 0.0 ?\n");
                        pamx[i][j] = 0.00001; /* by J. Thompson */
                    }
                    pamx[i][j] = log10(pamx[i][j]) * 1000.0;
                }

            if (disp) {
                reporterr(" after log\n");
                for (i = 0; i < 4; i++) {
                    for (j = 0; j < 4; j++)
                        reporterr("%+10.6f ", pamx[i][j]);
                    reporterr("\n");
                }
                reporterr("\n");
            }

            // ?????

            average = 0.0;
            for (i = 0; i < 4; i++)
                for (j = 0; j < 4; j++)
                    average += pamx[i][j] * freq[i] * freq[j];
            for (i = 0; i < 4; i++)
                for (j = 0; j < 4; j++)
                    pamx[i][j] -= average;

            average = 0.0;
            for (i = 0; i < 4; i++)
                average += pamx[i][i] * 1.0 / 4.0;

            for (i = 0; i < 4; i++)
                for (j = 0; j < 4; j++)
                    pamx[i][j] *= 600.0 / average;

            for (i = 0; i < 4; i++)
                for (j = 0; j < 4; j++)
                    pamx[i][j] -= offset;

            for (i = 0; i < 4; i++)
                for (j = 0; j < 4; j++)
                    pamx[i][j] = shishagonyuu(pamx[i][j]);

            if (disp) {
                reporterr(" after shishagonyuu\n");
                for (i = 0; i < 4; i++) {
                    for (j = 0; j < 4; j++)
                        reporterr("%+#6.10f", pamx[i][j]);
                    reporterr("\n");
                }
                reporterr("\n");
            }
            reporterr("done\n");
        }

        for (i = 0; i < 5; i++) {
            pamx[4][i] = pamx[3][i];
            pamx[i][4] = pamx[i][3];
        }

        for (i = 5; i < 10; i++)
            for (j = 5; j < 10; j++) {
                pamx[i][j] = pamx[i - 5][j - 5];
            }

        if (disp) {
            reporterr(" before dis\n");
            for (i = 0; i < 4; i++) {
                for (j = 0; j < 4; j++)
                    reporterr("%+#6.10f", pamx[i][j]);
                reporterr("\n");
            }
            reporterr("\n");
        }

        if (disp) {
            reporterr(" score matrix  \n");
            for (i = 0; i < 4; i++) {
                for (j = 0; j < 4; j++)
                    reporterr("%+#6.10f", pamx[i][j]);
                reporterr("\n");
            }
            reporterr("\n");
            exit(1);
        }

        for (i = 0; i < 26; i++)
            amino[i] = locaminon[i];
        for (i = 0; i < 26; i++)
            amino_grp[(int)amino[i]] = locgrpn[i];
        for (i = 0; i < 26; i++)
            for (j = 0; j < 26; j++)
                n_dis[i][j] = 0;
        for (i = 0; i < 10; i++)
            for (j = 0; j < 10; j++)
                n_dis[i][j] = shishagonyuu(pamx[i][j]);

        ambiguousscore(amino_n, n_dis);
        if (nwildcard)
            nscore(amino_n, n_dis);

        if (disp) {
            reporterr(" score matrix  \n");
            for (i = 0; i < 26; i++) {
                for (j = 0; j < 26; j++)
                    reporterr("%+6d", n_dis[i][j]);
                reporterr("\n");
            }
            reporterr("\n");
            reporterr("penalty = %d, penalty_ex = %d\n", penalty, penalty_ex);
            //exit( 1 );
        }

// RIBOSUM
#if 1
        average = 0.0;
        for (i = 0; i < 4; i++)
            for (j = 0; j < 4; j++)
                average += ribosum4[i][j] * freq[i] * freq[j];
        for (i = 0; i < 4; i++)
            for (j = 0; j < 4; j++)
                ribosum4[i][j] -= average;

        average = 0.0;
        for (i = 0; i < 4; i++)
            for (j = 0; j < 4; j++)
                for (k = 0; k < 4; k++)
                    for (m = 0; m < 4; m++) {
                        //			if( i%4==0&&j%4==3 || i%4==3&&j%4==0 || i%4==1&&j%4==2 || i%4==2&&j%4==1 || i%4==1&&j%4==3 || i%4==3&&j%4==1 )
                        //			if( k%4==0&&m%4==3 || k%4==3&&m%4==0 || k%4==1&&m%4==2 || k%4==2&&m%4==1 || k%4==1&&m%4==3 || k%4==3&&m%4==1 )
                        average += ribosum16[i * 4 + j][k * 4 + m] * freq[i] * freq[j] * freq[k] * freq[m];
                    }
        for (i = 0; i < 16; i++)
            for (j = 0; j < 16; j++)
                ribosum16[i][j] -= average;

        average = 0.0;
        for (i = 0; i < 4; i++)
            average += ribosum4[i][i] * freq[i];
        for (i = 0; i < 4; i++)
            for (j = 0; j < 4; j++)
                ribosum4[i][j] *= 600.0 / average;

        average = 0.0;
        average += ribosum16[0 * 4 + 3][0 * 4 + 3] * freq[0] * freq[3];  // AU
        average += ribosum16[3 * 4 + 0][3 * 4 + 0] * freq[3] * freq[0];  // UA
        average += ribosum16[1 * 4 + 2][1 * 4 + 2] * freq[1] * freq[2];  // CG
        average += ribosum16[2 * 4 + 1][2 * 4 + 1] * freq[2] * freq[1];  // GC
        average += ribosum16[1 * 4 + 3][1 * 4 + 3] * freq[1] * freq[3];  // GU
        average += ribosum16[3 * 4 + 1][3 * 4 + 1] * freq[3] * freq[1];  // UG
        for (i = 0; i < 16; i++)
            for (j = 0; j < 16; j++)
                ribosum16[i][j] *= 600.0 / average;

#if 1
        for (i = 0; i < 4; i++)
            for (j = 0; j < 4; j++)
                ribosum4[i][j] -= offset; /* extending gap cost ?????*/
        for (i = 0; i < 16; i++)
            for (j = 0; j < 16; j++)
                ribosum16[i][j] -= offset; /* extending gap cost ?????*/
#endif

        for (i = 0; i < 4; i++)
            for (j = 0; j < 4; j++)
                ribosum4[i][j] = shishagonyuu(ribosum4[i][j]);
        for (i = 0; i < 16; i++)
            for (j = 0; j < 16; j++)
                ribosum16[i][j] = shishagonyuu(ribosum16[i][j]);

        if (disp) {
            reporterr("ribosum after shishagonyuu\n");
            for (i = 0; i < 4; i++) {
                for (j = 0; j < 4; j++)
                    reporterr("%+#6.10f", ribosum4[i][j]);
                reporterr("\n");
            }
            reporterr("\n");
            reporterr("ribosum16 after shishagonyuu\n");
            for (i = 0; i < 16; i++) {
                for (j = 0; j < 16; j++)
                    reporterr("%+#7.0f", ribosum16[i][j]);
                reporterr("\n");
            }
            reporterr("\n");
        }
        //		reporterr(       "done\n" );

#if 1
        for (i = 0; i < 37; i++)
            for (j = 0; j < 37; j++)
                ribosumdis[i][j] = 0.0;  //iru
        for (m = 0; m < 9; m++)
            for (i = 0; i < 4; i++)  // loop
                for (k = 0; k < 9; k++)
                    for (j = 0; j < 4; j++)
                        ribosumdis[m * 4 + i][k * 4 + j] = ribosum4[i][j];  // loop-loop
        //			for( k=0; k<9; k++ ) for( j=0; j<4; j++ ) ribosumdis[m*4+i][k*4+j] = n_dis[i][j]; // loop-loop

        for (i = 0; i < 16; i++)
            for (j = 0; j < 16; j++)
                ribosumdis[i + 4][j + 4] = ribosum16[i][j];  // stem5-stem5
        for (i = 0; i < 16; i++)
            for (j = 0; j < 16; j++)
                ribosumdis[i + 20][j + 20] = ribosum16[i][j];  // stem5-stem5
#else  // do not use ribosum
        for (i = 0; i < 37; i++)
            for (j = 0; j < 37; j++)
                ribosumdis[i][j] = 0.0;  //iru
        for (m = 0; m < 9; m++)
            for (i = 0; i < 4; i++)  // loop
                for (k = 0; k < 9; k++)
                    for (j = 0; j < 4; j++)
                        ribosumdis[m * 4 + i][k * 4 + j] = n_dis[i][j];  // loop-loop
#endif

        if (disp) {
            reporterr("ribosumdis\n");
            for (i = 0; i < 37; i++) {
                for (j = 0; j < 37; j++)
                    reporterr("%+5d", ribosumdis[i][j]);
                reporterr("\n");
            }
            reporterr("\n");
        }
//		reporterr(       "done\n" );
#endif

        FreeDoubleMtx(pam1);
        FreeDoubleMtx(pamx);
        free(freq);

    } else if (dorp == 'p' && scoremtx == 1 && nblosum == -2) /* extended */
    {
        double* freq;
        double* freq1;
        double* datafreq;
        double  average;
        //		double tmp;
        double** n_distmp;
        int      userdefined;

        nalphabets = 0x100;
        nscoredalphabets = 0x100;
        charsize = 0x100;

        reporterr("nalphabets = %d\n", nalphabets);

        n_dis = AllocateIntMtx(nalphabets, nalphabets);
        n_disLN = AllocateDoubleMtx(nalphabets, nalphabets);
        n_distmp = AllocateDoubleMtx(nalphabets, nalphabets);
        datafreq = AllocateDoubleVec(nalphabets);
        freq = AllocateDoubleVec(nalphabets);

        if (ppenalty == NOTSPECIFIED)
            ppenalty = DEFAULTGOP_B;
        if (ppenalty_dist == NOTSPECIFIED)
            ppenalty_dist = ppenalty;
        if (ppenalty_OP == NOTSPECIFIED)
            ppenalty_OP = DEFAULTGOP_B;
        if (ppenalty_ex == NOTSPECIFIED)
            ppenalty_ex = DEFAULTGEP_B;
        if (ppenalty_EX == NOTSPECIFIED)
            ppenalty_EX = DEFAULTGEP_B;
        if (poffset == NOTSPECIFIED)
            poffset = DEFAULTOFS_B;
        if (pamN == NOTSPECIFIED)
            pamN = 0;
        if (kimuraR == NOTSPECIFIED)
            kimuraR = 1;
        penalty = (int)(600.0 / 1000.0 * ppenalty + 0.5);
        penalty_dist = (int)(600.0 / 1000.0 * ppenalty_dist + 0.5);
        penalty_shift = (int)(penalty_shift_factor * penalty);
        penalty_OP = (int)(600.0 / 1000.0 * ppenalty_OP + 0.5);
        penalty_ex = (int)(600.0 / 1000.0 * ppenalty_ex + 0.5);
        penalty_EX = (int)(600.0 / 1000.0 * ppenalty_EX + 0.5);
        offset = (int)(600.0 / 1000.0 * poffset + 0.5);
        offsetFFT = (int)(600.0 / 1000.0 * (-0) + 0.5);
        offsetLN = (int)(600.0 / 1000.0 * 100 + 0.5);
        penaltyLN = (int)(600.0 / 1000.0 * -2000 + 0.5);
        penalty_exLN = (int)(600.0 / 1000.0 * -100 + 0.5);

        userdefined = extendedmtx(n_distmp, freq, amino, amino_grp);

        if (trywarp)
            sprintf(shiftmodel, "%4.2f", -(double)penalty_shift / 600);
        else
            sprintf(shiftmodel, "noshift");

        sprintf(modelname, "Extended, %4.2f, %+4.2f, %+4.2f, %s", -(double)ppenalty / 1000, -(double)poffset / 1000, -(double)ppenalty_ex / 1000, shiftmodel);
#if 0
		for( i=0; i<26; i++ ) amino[i] = locaminod[i];
		for( i=0; i<26; i++ ) amino_grp[(int)amino[i]] = locgrpd[i];
		for( i=0; i<0x80; i++ ) amino_n[i] = 0;
		for( i=0; i<26; i++ ) amino_n[(int)amino[i]] = i;
#endif
        for (i = 0; i < 0x100; i++)
            amino_n[i] = -1;
        for (i = 0; i < nalphabets; i++) {
            amino_n[(unsigned char)amino[i]] = i;
            //			reporterr(       "i=%d, amino = %c, amino_n = %d\n", i, amino[i], amino_n[amino[i]] );
        }
        if (fmodel == 1) {
            calcfreq_extended(nseq, seq, datafreq);
            freq1 = datafreq;
        } else {
            calcfreq_from_scoremtx(n_distmp, datafreq);
            freq1 = datafreq;
        }
#if 1
        if (userdefined)
            if ((i = checkscoremtx(n_distmp, nseq, seq))) {
                reporterr("\n\nAlphabet %c (0x%x) is used in the sequence file but no score involving this alphabet is given in the matrix file.\n", i, i);
                reporterr("Check if the data is as intended.\n\n\n");
                exit(1);
            }
#endif

#if TEST
        reporterr("raw scoreing matrix : \n");
        for (i = 0; i < nalphabets; i++) {
            for (j = 0; j < nalphabets; j++) {
                fprintf(stdout, "%6.2f", n_distmp[i][j]);
            }
            fprintf(stdout, "\n");
        }
#endif
        if (fmodel == -1)
            average = 0.0;
        else {
#if TEST
            for (i = 0; i < nalphabets; i++)
                fprintf(stdout, "freq[%c] = %f, datafreq[%c] = %f, freq1[] = %f\n", amino[i], freq[i], amino[i], datafreq[i], freq1[i]);
#endif
            average = 0.0;
            for (i = 0; i < nalphabets; i++)
                for (j = 0; j < nalphabets; j++)
                    average += n_distmp[i][j] * freq1[i] * freq1[j];
        }
#if 0
		if( disp ) fprintf( stdout, "####### average2  = %f\n", average );
#endif

        for (i = 0; i < nalphabets; i++)
            for (j = 0; j < nalphabets; j++) {
                if (n_distmp[i][j] == -1.0)
                    n_distmp[i][j] = -average;
                else
                    n_distmp[i][j] -= average;
            }
#if 0
		fprintf( stderr, "average2 = %f\n", average );
		fprintf( stderr, "after average subtruction : \n" );
		for( i=0; i<nalphabets; i++ )
		{
			fprintf( stderr, "i=%d, %x\n", i, i );
			for( j=0; j<nalphabets; j++ ) 
			{
				fprintf( stderr, "%6.2f", n_distmp[i][j] );
			}
			fprintf( stderr, "\n" );
		}
#endif

        average = 0.0;
        for (i = 0; i < nalphabets; i++)
            average += n_distmp[i][i] * freq1[i];
#if 0
		if( disp ) fprintf( stdout, "####### average1  = %f\n", average );
#endif

        if (average < 0.0) {
            reporterr("\nUnrealistic scoring matrix.  Give larger positive values to matches (A/A, B/B, etc).\n\n");
            exit(1);
        }

        for (i = 0; i < nalphabets; i++)
            for (j = 0; j < nalphabets; j++)
                n_distmp[i][j] *= 600.0 / average;
#if TEST
        fprintf(stderr, "after average division : \n");
        for (i = 0; i < nalphabets; i++) {
            fprintf(stderr, "i=%d, %x\n", i, i);
            for (j = 0; j <= i; j++) {
                fprintf(stderr, "%7.1f", n_distmp[i][j]);
            }
            fprintf(stderr, "\n");
        }
#endif

        for (i = 0; i < nalphabets; i++)
            for (j = 0; j < nalphabets; j++)
                n_distmp[i][j] -= offset;
#if TEST
        fprintf(stderr, "after offset subtruction (offset = %d): \n", offset);
        for (i = 0; i < nalphabets; i++) {
            fprintf(stderr, "i=%d, %x\n", i, i);
            for (j = 0; j <= i; j++) {
                fprintf(stderr, "%30.10f", n_distmp[i][j]);
            }
            fprintf(stderr, "\n");
        }
#endif
#if 0
/* ���� �������������������� */
			penalty -= offset;
#endif

        for (i = 0; i < nalphabets; i++)
            for (j = 0; j < nalphabets; j++)
                n_distmp[i][j] = shishagonyuu(n_distmp[i][j]);

        if (disp) {
            fprintf(stdout, "freq = \n");
            for (i = 0; i < nalphabets; i++)
                fprintf(stdout, "%c %f\n", amino[i], freq1[i]);
            fprintf(stdout, " scoring matrix  \n");
            for (i = 0; i < nalphabets; i++) {
                fprintf(stdout, "%c    ", amino[i]);
                for (j = 0; j < nalphabets; j++)
                    fprintf(stdout, "%5.0f", n_distmp[i][j]);
                fprintf(stdout, "\n");
            }
            fprintf(stdout, "     ");
            for (i = 0; i < nalphabets; i++)
                fprintf(stdout, "    %c", amino[i]);

            average = 0.0;
            for (i = 0; i < nalphabets; i++)
                for (j = 0; j < nalphabets; j++)
                    average += n_distmp[i][j] * freq1[i] * freq1[j];
            fprintf(stdout, "average = %f\n", average);

            average = 0.0;
            for (i = 0; i < nalphabets; i++)
                average += n_distmp[i][i] * freq1[i];
            fprintf(stdout, "itch average = %f\n", average);
            reporterr("parameters: %d, %d, %d\n", penalty, penalty_ex, offset);

            exit(1);
        }

        for (i = 0; i < nalphabets; i++)
            for (j = 0; j < nalphabets; j++)
                n_dis[i][j] = 0;
        for (i = 0; i < nalphabets; i++)
            for (j = 0; j < nalphabets; j++)
                n_dis[i][j] = (int)n_distmp[i][j];
        for (i = 0; i < nalphabets; i++)
            for (j = 0; j < nalphabets; j++)
                n_dis[i][amino_n['-']] = n_dis[amino_n['-']][i] = 0;

        FreeDoubleMtx(n_distmp);
        FreeDoubleVec(datafreq);
        FreeDoubleVec(freq);

        //		reporterr(       "done.\n" );

    } else if (dorp == 'p' && scoremtx == 1) /* Blosum, user-defined */
    {
        double* freq;
        double* freq1;
        double* datafreq;
        double  average;
        double  iaverage;
        //		double tmp;
        double** n_distmp;
        int      rescale = 1;

        if (nblosum == 0) {
            reporterr("nblosum=%d??\n", nblosum);
            exit(1);
        }
        //		if( nblosum < 0 )
        //		{
        //			nblosum *= -1;
        //			makeaverage0 = 0;
        //		}

        nalphabets = 26;
        nscoredalphabets = 20;
        charsize = 0x80;

        n_dis = AllocateIntMtx(nalphabets, nalphabets);
        n_disLN = AllocateDoubleMtx(nalphabets, nalphabets);
        n_distmp = AllocateDoubleMtx(20, 20);
        datafreq = AllocateDoubleVec(20);
        freq = AllocateDoubleVec(20);

        if (ppenalty == NOTSPECIFIED)
            ppenalty = DEFAULTGOP_B;
        if (ppenalty_dist == NOTSPECIFIED)
            ppenalty_dist = ppenalty;
        if (ppenalty_OP == NOTSPECIFIED)
            ppenalty_OP = DEFAULTGOP_B;
        if (ppenalty_ex == NOTSPECIFIED)
            ppenalty_ex = DEFAULTGEP_B;
        if (ppenalty_EX == NOTSPECIFIED)
            ppenalty_EX = DEFAULTGEP_B;
        if (poffset == NOTSPECIFIED)
            poffset = DEFAULTOFS_B;
        if (pamN == NOTSPECIFIED)
            pamN = 0;
        if (kimuraR == NOTSPECIFIED)
            kimuraR = 1;
        penalty = (int)(600.0 / 1000.0 * ppenalty + 0.5);
        penalty_dist = (int)(600.0 / 1000.0 * ppenalty_dist + 0.5);
        penalty_shift = (int)(penalty_shift_factor * penalty);
        penalty_OP = (int)(600.0 / 1000.0 * ppenalty_OP + 0.5);
        penalty_ex = (int)(600.0 / 1000.0 * ppenalty_ex + 0.5);
        penalty_EX = (int)(600.0 / 1000.0 * ppenalty_EX + 0.5);
        offset = (int)(600.0 / 1000.0 * poffset + 0.5);
        offsetFFT = (int)(600.0 / 1000.0 * (-0) + 0.5);
        offsetLN = (int)(600.0 / 1000.0 * 100 + 0.5);
        penaltyLN = (int)(600.0 / 1000.0 * -2000 + 0.5);
        penalty_exLN = (int)(600.0 / 1000.0 * -100 + 0.5);

        BLOSUMmtx(nblosum, n_distmp, freq, amino, amino_grp, &rescale);

        reporterr("rescale = %d\n", rescale);

        if (trywarp)
            sprintf(shiftmodel, "%4.2f", -(double)penalty_shift / 600);
        else
            sprintf(shiftmodel, "noshift");

        if (nblosum == -1)
            sprintf(modelname, "User-defined, %4.2f, %+4.2f, %+4.2f, %s", -(double)ppenalty / 1000, -(double)poffset / 1000, -(double)ppenalty_ex / 1000, shiftmodel);
        else
            sprintf(modelname, "BLOSUM%d, %4.2f, %+4.2f, %+4.2f, %s", nblosum, -(double)ppenalty / 1000, -(double)poffset / 1000, -(double)ppenalty_ex / 1000, shiftmodel);
#if 0
		for( i=0; i<26; i++ ) amino[i] = locaminod[i];
		for( i=0; i<26; i++ ) amino_grp[(int)amino[i]] = locgrpd[i];
		for( i=0; i<0x80; i++ ) amino_n[i] = 0;
		for( i=0; i<26; i++ ) amino_n[(int)amino[i]] = i;
#endif
        for (i = 0; i < 0x80; i++)
            amino_n[i] = -1;
        for (i = 0; i < 26; i++)
            amino_n[(int)amino[i]] = i;
        if (fmodel == 1) {
            calcfreq(nseq, seq, datafreq);
            freq1 = datafreq;
        } else
            freq1 = freq;
#if TEST
        reporterr("raw scoreing matrix : \n");
        for (i = 0; i < 20; i++) {
            for (j = 0; j < 20; j++) {
                fprintf(stdout, "%6.2f", n_distmp[i][j]);
            }
            fprintf(stdout, "\n");
        }
#endif
        if (fmodel == -1)
            average = 0.0;
        else {
#if TEST
            for (i = 0; i < 20; i++)
                fprintf(stdout, "freq[%c] = %f, datafreq[%c] = %f, freq1[] = %f\n", amino[i], freq[i], amino[i], datafreq[i], freq1[i]);
#endif
            average = 0.0;
            for (i = 0; i < 20; i++)
                for (j = 0; j < 20; j++)
                    average += n_distmp[i][j] * freq1[i] * freq1[j];
        }
#if TEST
        if (disp)
            fprintf(stdout, "####### average2  = %f\n", average);
#endif

        if (rescale) {
            for (i = 0; i < 20; i++)
                for (j = 0; j < 20; j++)
                    n_distmp[i][j] -= average;
        }
#if TEST
        fprintf(stdout, "average2 = %f\n", average);
        fprintf(stdout, "after average subtraction : \n");
        for (i = 0; i < 20; i++) {
            for (j = 0; j < 20; j++) {
                fprintf(stdout, "%6.2f", n_distmp[i][j]);
            }
            fprintf(stdout, "\n");
        }
#endif

        average = 0.0;
        for (i = 0; i < 20; i++)
            average += n_distmp[i][i] * freq1[i];
#if TEST
        if (disp)
            fprintf(stdout, "####### average1  = %f\n", average);
#endif

        if (rescale) {
            for (i = 0; i < 20; i++)
                for (j = 0; j < 20; j++)
                    n_distmp[i][j] *= 600.0 / average;
        } else {
            for (i = 0; i < 20; i++)
                for (j = 0; j < 20; j++)
                    n_distmp[i][j] *= 600.0;
        }
#if TEST
        fprintf(stdout, "after average division : \n");
        for (i = 0; i < 20; i++) {
            for (j = 0; j <= i; j++) {
                fprintf(stdout, "%7.1f", n_distmp[i][j]);
            }
            fprintf(stdout, "\n");
        }
#endif

        for (i = 0; i < 20; i++)
            for (j = 0; j < 20; j++)
                n_distmp[i][j] -= offset;
#if TEST
        fprintf(stdout, "after offset substruction (offset = %d): \n", offset);
        for (i = 0; i < 20; i++) {
            for (j = 0; j <= i; j++) {
                fprintf(stdout, "%7.1f", n_distmp[i][j]);
            }
            fprintf(stdout, "\n");
        }
#endif
#if 0
/* ���� �������������������� */
			penalty -= offset;
#endif

        for (i = 0; i < 20; i++)
            for (j = 0; j < 20; j++)
                n_distmp[i][j] = shishagonyuu(n_distmp[i][j]);

        if (disp) {
            fprintf(stderr, " scoring matrix  \n");
            for (i = 0; i < 20; i++) {
                fprintf(stderr, "%c    ", amino[i]);
                for (j = 0; j < 20; j++)
                    fprintf(stderr, "%5.0f", n_distmp[i][j]);
                fprintf(stderr, "\n");
            }
            fprintf(stderr, "     ");
            for (i = 0; i < 20; i++)
                fprintf(stderr, "    %c", amino[i]);

            average = 0.0;
            for (i = 0; i < 20; i++)
                for (j = 0; j < 20; j++)
                    average += n_distmp[i][j] * freq1[i] * freq1[j];
            fprintf(stderr, "\naverage = %f\n", average);

            iaverage = 0.0;
            for (i = 0; i < 20; i++)
                iaverage += n_distmp[i][i] * freq1[i];
            fprintf(stderr, "itch average = %f, E=%f\n", iaverage, average / iaverage);
            reporterr("parameters: %d, %d, %d\n", penalty, penalty_ex, offset);

            exit(1);
        }

        for (i = 0; i < 26; i++)
            for (j = 0; j < 26; j++)
                n_dis[i][j] = 0;
        for (i = 0; i < 20; i++)
            for (j = 0; j < 20; j++)
                n_dis[i][j] = (int)n_distmp[i][j];

        FreeDoubleMtx(n_distmp);
        FreeDoubleVec(datafreq);
        FreeDoubleVec(freq);

        //		reporterr(       "done.\n" );

    } else if (dorp == 'p' && scoremtx == 2) /* Miyata-Yasunaga */
    {
        reporterr("Not supported\n");
        exit(1);
    } else /* JTT */
    {
        double** rsr;
        double** pam1;
        double** pamx;
        double*  freq;
        double*  freq1;
        double*  mutab;
        double*  datafreq;
        double   average;
        double   iaverage;
        double   tmp;
        double   delta;
        int      makeaverage0;

        nalphabets = 26;
        nscoredalphabets = 20;
        charsize = 0x80;

        n_dis = AllocateIntMtx(nalphabets, nalphabets);
        n_disLN = AllocateDoubleMtx(nalphabets, nalphabets);
        rsr = AllocateDoubleMtx(20, 20);
        pam1 = AllocateDoubleMtx(20, 20);
        pamx = AllocateDoubleMtx(20, 20);
        freq = AllocateDoubleVec(20);
        mutab = AllocateDoubleVec(20);
        datafreq = AllocateDoubleVec(20);

        if (ppenalty == NOTSPECIFIED)
            ppenalty = DEFAULTGOP_J;
        if (ppenalty_dist == NOTSPECIFIED)
            ppenalty_dist = ppenalty;
        if (ppenalty_OP == NOTSPECIFIED)
            ppenalty_OP = DEFAULTGOP_J;
        if (ppenalty_ex == NOTSPECIFIED)
            ppenalty_ex = DEFAULTGEP_J;
        if (ppenalty_EX == NOTSPECIFIED)
            ppenalty_EX = DEFAULTGEP_J;
        if (poffset == NOTSPECIFIED)
            poffset = DEFAULTOFS_J;
        if (pamN == NOTSPECIFIED)
            pamN = DEFAULTPAMN;
        if (kimuraR == NOTSPECIFIED)
            kimuraR = 1;

        if (pamN == 0) {
            reporterr("pamN=%d??\n", pamN);
            exit(1);
        }
        if (pamN < 0) {
            pamN *= -1;
            makeaverage0 = 0;
        } else {
            makeaverage0 = 1;
        }

        penalty = (int)(600.0 / 1000.0 * ppenalty + 0.5);
        penalty_dist = (int)(600.0 / 1000.0 * ppenalty_dist + 0.5);
        penalty_shift = (int)(penalty_shift_factor * penalty);
        penalty_OP = (int)(600.0 / 1000.0 * ppenalty_OP + 0.5);
        penalty_ex = (int)(600.0 / 1000.0 * ppenalty_ex + 0.5);
        penalty_EX = (int)(600.0 / 1000.0 * ppenalty_EX + 0.5);
        offset = (int)(600.0 / 1000.0 * poffset + 0.5);
        offsetFFT = (int)(600.0 / 1000.0 * (-0) + 0.5);
        offsetLN = (int)(600.0 / 1000.0 * 100 + 0.5);
        penaltyLN = (int)(600.0 / 1000.0 * -2000 + 0.5);
        penalty_exLN = (int)(600.0 / 1000.0 * -100 + 0.5);

        if (trywarp)
            sprintf(shiftmodel, "%4.2f", -(double)penalty_shift / 600);
        else
            sprintf(shiftmodel, "noshift");

        sprintf(modelname, "%s %dPAM, %4.2f, %4.2f, %s", (TMorJTT == TM) ? "Transmembrane" : "JTT", pamN, -(double)ppenalty / 1000, -(double)poffset / 1000, shiftmodel);

        JTTmtx(rsr, freq, amino, amino_grp, (int)(TMorJTT == TM));

        for (i = 0; i < 0x80; i++)
            amino_n[i] = -1;
        for (i = 0; i < 26; i++)
            amino_n[(int)amino[i]] = i;
        if (fmodel == 1) {
            calcfreq(nseq, seq, datafreq);
            freq1 = datafreq;
        } else
            freq1 = freq;

#if TEST
        fprintf(stdout, "rsr = \n");
        for (i = 0; i < 20; i++) {
            for (j = 0; j < 20; j++) {
                fprintf(stdout, "%9.2f ", rsr[i][j]);
            }
            fprintf(stdout, "\n");
        }
#endif

        reporterr("generating %dPAM %s scoring matrix for amino acids ... ", pamN, (TMorJTT == TM) ? "Transmembrane" : "JTT");

        tmp = 0.0;
        for (i = 0; i < 20; i++) {
            mutab[i] = 0.0;
            for (j = 0; j < 20; j++)
                mutab[i] += rsr[i][j] * freq1[j];
            tmp += mutab[i] * freq1[i];
        }
#if TEST
        fprintf(stdout, "mutability = \n");
        for (i = 0; i < 20; i++)
            fprintf(stdout, "%5.3f\n", mutab[i]);

        fprintf(stdout, "tmp = %f\n", tmp);
#endif
        delta = 0.01 / tmp;
        for (i = 0; i < 20; i++) {
            for (j = 0; j < 20; j++) {
                if (i != j)
                    pam1[i][j] = delta * rsr[i][j] * freq1[j];
                else
                    pam1[i][j] = 1.0 - delta * mutab[i];
            }
        }

        if (disp) {
            fprintf(stdout, "pam1 = \n");
            for (i = 0; i < 20; i++) {
                for (j = 0; j < 20; j++) {
                    fprintf(stdout, "%9.6f ", pam1[i][j]);
                }
                fprintf(stdout, "\n");
            }
        }

        MtxuntDouble(pamx, 20);
        for (x = 0; x < pamN; x++)
            MtxmltDouble(pamx, pam1, 20);

        for (i = 0; i < 20; i++)
            for (j = 0; j < 20; j++)
                pamx[i][j] /= freq1[j];

        for (i = 0; i < 20; i++)
            for (j = 0; j < 20; j++) {
                if (pamx[i][j] == 0.0) {
                    reporterr("WARNING: pamx[%d][%d] = 0.0?\n", i, j);
                    pamx[i][j] = 0.00001; /* by J. Thompson */
                }
                pamx[i][j] = log10(pamx[i][j]) * 1000.0;
            }

#if TEST
        fprintf(stdout, "raw scoring matrix : \n");
        for (i = 0; i < 20; i++) {
            for (j = 0; j < 20; j++) {
                fprintf(stdout, "%5.0f", pamx[i][j]);
            }
            fprintf(stdout, "\n");
        }
        average = tmp = 0.0;
        for (i = 0; i < 20; i++)
            for (j = 0; j < 20; j++) {
                average += pamx[i][j] * freq1[i] * freq1[j];
                tmp += freq1[i] * freq1[j];
            }
        average /= tmp;
        fprintf(stdout, "Zenbu average = %f, tmp = %f \n", average, tmp);
        average = tmp = 0.0;
        for (i = 0; i < 20; i++)
            for (j = i; j < 20; j++) {
                average += pamx[i][j] * freq1[i] * freq1[j];
                tmp += freq1[i] * freq1[j];
            }
        average /= tmp;
        fprintf(stdout, "Zenbu average2 = %f, tmp = %f \n", average, tmp);
        average = tmp = 0.0;
        for (i = 0; i < 20; i++) {
            average += pamx[i][i] * freq1[i];
            tmp += freq1[i];
        }
        average /= tmp;
        fprintf(stdout, "Itch average = %f, tmp = %f \n", average, tmp);
#endif

#if NORMALIZE1
        if (fmodel == -1)
            average = 0.0;
        else {
#if TEST
            for (i = 0; i < 20; i++)
                fprintf(stdout, "freq[%c] = %f, datafreq[%c] = %f, freq1[] = %f\n", amino[i], freq[i], amino[i], datafreq[i], freq1[i]);
#endif
            average = 0.0;
            for (i = 0; i < 20; i++)
                for (j = 0; j < 20; j++)
                    average += pamx[i][j] * freq1[i] * freq1[j];
        }
#if TEST
        fprintf(stdout, "####### average2  = %f\n", average);
#endif

        if (makeaverage0) {
            for (i = 0; i < 20; i++)
                for (j = 0; j < 20; j++)
                    pamx[i][j] -= average;
        }
#if TEST
        fprintf(stdout, "average2 = %f\n", average);
        fprintf(stdout, "after average substruction : \n");
        for (i = 0; i < 20; i++) {
            for (j = 0; j < 20; j++) {
                fprintf(stdout, "%5.0f", pamx[i][j]);
            }
            fprintf(stdout, "\n");
        }
#endif

        average = 0.0;
        for (i = 0; i < 20; i++)
            average += pamx[i][i] * freq1[i];
#if TEST
        fprintf(stdout, "####### average1  = %f\n", average);
#endif

        for (i = 0; i < 20; i++)
            for (j = 0; j < 20; j++)
                pamx[i][j] *= 600.0 / average;
#if TEST
        fprintf(stdout, "after average division : \n");
        for (i = 0; i < 20; i++) {
            for (j = 0; j <= i; j++) {
                fprintf(stdout, "%5.0f", pamx[i][j]);
            }
            fprintf(stdout, "\n");
        }
#endif

        for (i = 0; i < 20; i++)
            for (j = 0; j < 20; j++)
                pamx[i][j] -= offset;
#if TEST
        fprintf(stdout, "after offset substruction (offset = %d): \n", offset);
        for (i = 0; i < 20; i++) {
            for (j = 0; j <= i; j++) {
                fprintf(stdout, "%5.0f", pamx[i][j]);
            }
            fprintf(stdout, "\n");
        }
#endif
#if 0
/* ���� �������������������� */
			penalty -= offset;
#endif

        for (i = 0; i < 20; i++)
            for (j = 0; j < 20; j++)
                pamx[i][j] = shishagonyuu(pamx[i][j]);

#else

        average = 0.0;
        for (i = 0; i < 20; i++)
            for (j = 0; j < 20; j++)
                average += pamx[i][j];
        average /= 400.0;

        for (i = 0; i < 20; i++)
            for (j = 0; j < 20; j++) {
                pamx[i][j] -= average;
                pamx[i][j] = shishagonyuu(pamx[i][j]);
            }
#endif
        if (disp) {
            fprintf(stdout, " scoring matrix  \n");
            for (i = 0; i < 20; i++) {
                fprintf(stdout, "%c    ", amino[i]);
                for (j = 0; j < 20; j++)
                    fprintf(stdout, "%5.0f", pamx[i][j]);
                fprintf(stdout, "\n");
            }
            fprintf(stdout, "     ");
            for (i = 0; i < 20; i++)
                fprintf(stdout, "    %c", amino[i]);

            average = 0.0;
            for (i = 0; i < 20; i++)
                for (j = 0; j < 20; j++)
                    average += pamx[i][j] * freq1[i] * freq1[j];
            fprintf(stdout, "\naverage = %f\n", average);

            iaverage = 0.0;
            for (i = 0; i < 20; i++)
                iaverage += pamx[i][i] * freq1[i];
            fprintf(stdout, "itch average = %f, E=%f\n", average, average / iaverage);
            reporterr("parameters: %d, %d, %d\n", penalty, penalty_ex, offset);

            exit(1);
        }

        for (i = 0; i < 26; i++)
            for (j = 0; j < 26; j++)
                n_dis[i][j] = 0;
        for (i = 0; i < 20; i++)
            for (j = 0; j < 20; j++)
                n_dis[i][j] = (int)pamx[i][j];

        reporterr("done.\n");
        FreeDoubleMtx(rsr);
        FreeDoubleMtx(pam1);
        FreeDoubleMtx(pamx);
        FreeDoubleVec(freq);
        FreeDoubleVec(mutab);
        FreeDoubleVec(datafreq);
    }

#if DEBUG
    reporterr("scoremtx = %d\n", scoremtx);
    reporterr("amino[] = %s\n", amino);
#endif

    amino_dis = AllocateIntMtx(charsize, charsize);
    amino_dis_consweight_multi = AllocateDoubleMtx(charsize, charsize);

    //	reporterr( "charsize=%d\n", charsize );

    for (i = 0; i < charsize; i++)
        amino_n[i] = -1;
    for (i = 0; i < nalphabets; i++)
        amino_n[(int)amino[i]] = i;
    for (i = 0; i < charsize; i++)
        for (j = 0; j < charsize; j++)
            amino_dis[i][j] = 0;
    for (i = 0; i < nalphabets; i++)
        for (j = 0; j < nalphabets; j++)
            n_disLN[i][j] = 0;
    for (i = 0; i < charsize; i++)
        for (j = 0; j < charsize; j++)
            amino_dis_consweight_multi[i][j] = 0.0;

    n_dis_consweight_multi = AllocateDoubleMtx(nalphabets, nalphabets);
    n_disFFT = AllocateIntMtx(nalphabets, nalphabets);
    for (i = 0; i < nalphabets; i++)
        for (j = 0; j < nalphabets; j++) {
            amino_dis[(int)amino[i]][(int)amino[j]] = n_dis[i][j];
            n_dis_consweight_multi[i][j] = (double)n_dis[i][j] * consweight_multi;
            amino_dis_consweight_multi[(int)amino[i]][(int)amino[j]] = (double)n_dis[i][j] * consweight_multi;
        }

    if (dorp == 'd') /* DNA */
    {
#if 0  // ???
	    for( i=0; i<5; i++) for( j=0; j<5; j++ )
        	n_disLN[i][j] = (double)n_dis[i][j] + offset - offsetLN;
	    for( i=5; i<10; i++) for( j=5; j<10; j++ )
        	n_disLN[i][j] = (double)n_dis[i][j] + offset - offsetLN;
	    for( i=0; i<5; i++) for( j=0; j<5; j++ )
        	n_disFFT[i][j] = n_dis[i][j] + offset - offsetFFT;
	    for( i=5; i<10; i++) for( j=5; j<10; j++ )
        	n_disFFT[i][j] = n_dis[i][j] + offset - offsetFFT;
#else
        for (i = 0; i < 10; i++)
            for (j = 0; j < 10; j++)
                n_disLN[i][j] = (double)n_dis[i][j] + offset - offsetLN;
        for (i = 0; i < 10; i++)
            for (j = 0; j < 10; j++)
                n_disFFT[i][j] = n_dis[i][j] + offset - offsetFFT;
#endif
    } else  // protein
    {
        for (i = 0; i < 20; i++)
            for (j = 0; j < 20; j++)
                n_disLN[i][j] = (double)n_dis[i][j] + offset - offsetLN;
        for (i = 0; i < 20; i++)
            for (j = 0; j < 20; j++)
                n_disFFT[i][j] = n_dis[i][j] + offset - offsetFFT;
    }

#if 0
		reporterr(       "amino_dis (offset = %d): \n", offset );
		for( i=0; i<20; i++ )
		{
			for( j=0; j<20; j++ ) 
			{
				reporterr(       "%5d", amino_dis[(int)amino[i]][(int)amino[j]] );
			}
			reporterr(       "\n" );
		}

		reporterr(       "amino_disLN (offsetLN = %d): \n", offsetLN );
		for( i=0; i<20; i++ )
		{
			for( j=0; j<20; j++ ) 
			{
				reporterr(       "%5d", amino_disLN[(int)amino[i]][(int)amino[j]] );
			}
			reporterr(       "\n" );
		}

		reporterr(       "n_dis (offset = %d): \n", offset );
		for( i=0; i<26; i++ )
		{
			for( j=0; j<26; j++ ) 
			{
				reporterr(       "%5d", n_dis[i][j] );
			}
			reporterr(       "\n" );
		}

		reporterr(       "n_disFFT (offsetFFT = %d): \n", offsetFFT );
		for( i=0; i<26; i++ )
		{
			for( j=0; j<26; j++ ) 
			{
				reporterr(       "%5d", n_disFFT[i][j] );
			}
			reporterr(       "\n" );
		}
exit( 1 );
#endif

    ppid = 0;

    if (fftThreshold == NOTSPECIFIED) {
        fftThreshold = FFT_THRESHOLD;
    }
    if (fftWinSize == NOTSPECIFIED) {
        if (dorp == 'd')
            fftWinSize = FFT_WINSIZE_D;
        else
            fftWinSize = FFT_WINSIZE_P;
    }

    if (fftscore) {
        double av, sd;

        for (i = 0; i < 20; i++)
            polarity[i] = polarity_[i];
        for (av = 0.0, i = 0; i < 20; i++)
            av += polarity[i];
        av /= 20.0;
        for (sd = 0.0, i = 0; i < 20; i++)
            sd += (polarity[i] - av) * (polarity[i] - av);
        sd /= 20.0;
        sd = sqrt(sd);
        for (i = 0; i < 20; i++)
            polarity[i] -= av;
        for (i = 0; i < 20; i++)
            polarity[i] /= sd;

        for (i = 0; i < 20; i++)
            volume[i] = volume_[i];
        for (av = 0.0, i = 0; i < 20; i++)
            av += volume[i];
        av /= 20.0;
        for (sd = 0.0, i = 0; i < 20; i++)
            sd += (volume[i] - av) * (volume[i] - av);
        sd /= 20.0;
        sd = sqrt(sd);
        for (i = 0; i < 20; i++)
            volume[i] -= av;
        for (i = 0; i < 20; i++)
            volume[i] /= sd;

#if 0
		for( i=0; i<20; i++ ) fprintf( stdout, "amino=%c, pol = %f<-%f, vol = %f<-%f\n", amino[i], polarity[i], polarity_[i], volume[i], volume_[i] );
		for( i=0; i<20; i++ ) fprintf( stdout, "%c %+5.3f %+5.3f\n", amino[i], volume[i], polarity[i] );
#endif
    }
}

void
freeconstants() {
    if (n_disLN)
        FreeDoubleMtx(n_disLN);
    n_disLN = NULL;
    if (n_dis)
        FreeIntMtx(n_dis);
    n_dis = NULL;
    if (n_disFFT)
        FreeIntMtx(n_disFFT);
    n_disFFT = NULL;
    if (n_dis_consweight_multi)
        FreeDoubleMtx(n_dis_consweight_multi);
    n_dis_consweight_multi = NULL;
    if (amino_dis)
        FreeIntMtx(amino_dis);
    amino_dis = NULL;
    if (amino_dis_consweight_multi)
        FreeDoubleMtx(amino_dis_consweight_multi);
    amino_dis_consweight_multi = NULL;
}
