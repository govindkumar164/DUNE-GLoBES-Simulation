#include <stdio.h>
#include <math.h>
#include <globes/globes.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main(int argc, char *argv[])
{
    /* Initialize GLoBES */
    glbInit(argv[0]);
    glbInitExperiment("DUNE_GLoBES.glb",
                      &glb_experiment_list[0],
                      &glb_num_of_exps);

    /* Parameter containers */
    glb_params true_values = glbAllocParams();
    glb_params test_values = glbAllocParams();
    glb_params input_errors = glbAllocParams();
    glb_params minimum = glbAllocParams();

    /* Prior uncertainties */
    glbDefineParams(input_errors,
                     0.5878 * 0.0215,   /* theta12 */
                     0.1487 * 0.0135,   /* theta13 */
                     0.8465 * 0.0326,   /* theta23 */
                     0.0,               /* deltaCP */
                     7.49e-5 * 0.0251,  /* dm21 */
                     2.534e-3 * 0.0094);/* dm31 */

    glbSetDensityParams(input_errors, 0.05, GLB_ALL);

    /* theta13 values (radians) */
    double theta13_values[3] = {
        0.1428,   /* 8.18° */
        0.1487,   /* 8.52° best fit */
        0.1548    /* 8.87° */
    };

    FILE *fp = fopen("CPV_theta13.dat", "w");

    /* Loop over deltaCP */
    for (int i = 0; i <= 100; i++)
    {
        double deltacp = -M_PI + i * (2 * M_PI / 100.0);
        double sigma[3];

        /* Loop over theta13 */
        for (int j = 0; j < 3; j++)
        {
            double theta13 = theta13_values[j];

            /* TRUE parameters */
            glbDefineParams(true_values,
                             0.5878,
                             theta13,
                             0.8465,
                             deltacp,
                             7.49e-5,
                             2.534e-3);

            glbSetDensityParams(true_values, 1.0, GLB_ALL);

            /* Generate Asimov data */
            glbSetOscillationParameters(true_values);
            glbSetRates();
            glbSetCentralValues(true_values);
            glbSetInputErrors(input_errors);

            /* TRUE chi2 */
            glb_projection proj_true = glbAllocProjection();

            glbDefineProjection(proj_true,
                                 GLB_FREE, GLB_FREE, GLB_FREE,
                                 GLB_FREE, GLB_FREE, GLB_FREE);

            glbSetProjection(proj_true);

            double chi2_true =
                glbChiNP(true_values, minimum, GLB_ALL);

            glbFreeProjection(proj_true);

            /* Test delta = 0, pi */
            glb_projection proj_test =
                glbAllocProjection();

            glbDefineProjection(proj_test,
                                 GLB_FREE,
                                 GLB_FREE,
                                 GLB_FREE,
                                 GLB_FIXED,
                                 GLB_FREE,
                                 GLB_FREE);

            glbSetProjection(proj_test);
            glbCopyParams(true_values, test_values);

            /* delta = 0 */
            glbSetOscParams(test_values, 0.0, GLB_DELTA_CP);
            double chi2_0 =
                glbChiNP(test_values, minimum, GLB_ALL);

            /* delta = pi */
            glbSetOscParams(test_values, M_PI, GLB_DELTA_CP);
            double chi2_pi =
                glbChiNP(test_values, minimum, GLB_ALL);

            glbFreeProjection(proj_test);

            /* Sensitivity */
            double chi2_test = fmin(chi2_0, chi2_pi);
            double delta_chi2 = chi2_test - chi2_true;

            if (delta_chi2 < 0)
                delta_chi2 = 0;

            sigma[j] = sqrt(delta_chi2);
        }

        fprintf(fp, "%f %f %f %f\n",
                deltacp,
                sigma[0],
                sigma[1],
                sigma[2]);
    }

    fclose(fp);

    /* Free memory */
    glbFreeParams(true_values);
    glbFreeParams(test_values);
    glbFreeParams(input_errors);
    glbFreeParams(minimum);

    printf("Finished successfully.\n");
    return 0;
}