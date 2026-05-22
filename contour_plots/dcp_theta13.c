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

    /* NuFIT 6.0 priors */
    glbDefineParams(input_errors,
                     0.5878 * 0.0215,   /* theta12 */
                     0.1487 * 0.0135,   /* theta13 */
                     0.8465 * 0.0326,   /* theta23 */
                     0.0,               /* deltaCP */
                     7.49e-5 * 0.0251,  /* dm21 */
                     2.534e-3 * 0.0094);/* dm31 */

    glbSetDensityParams(input_errors, 0.05, GLB_ALL);

    /* True values (NuFIT 6.0) */
    double delta_true = -0.79 * M_PI;
    double sin2_2theta13_true = 0.0859;

    /* Convert sin²(2θ13) → θ13 */
    double theta13_true =
        0.5 * asin(sqrt(sin2_2theta13_true));

    glbDefineParams(true_values,
                     0.5878,
                     theta13_true,
                     0.8465,
                     delta_true,
                     7.49e-5,
                     2.534e-3);

    glbSetDensityParams(true_values, 1.0, GLB_ALL);

    /* Generate Asimov data */
    glbSetOscillationParameters(true_values);
    glbSetRates();
    glbSetCentralValues(true_values);
    glbSetInputErrors(input_errors);
    glbCopyParams(true_values, minimum);

    /* Global minimum χ² */
    glb_projection proj_min =
        glbAllocProjection();

    glbDefineProjection(proj_min,
                         GLB_FREE, GLB_FREE, GLB_FREE,
                         GLB_FREE, GLB_FREE, GLB_FREE);

    glbSetProjection(proj_min);

    double chi2_min =
        glbChiNP(true_values, minimum, GLB_ALL);

    glbFreeProjection(proj_min);

    /* Open output file */
    FILE *fp = fopen("contour_13yr.dat", "w");

    /* Scan grid */
    int n_delta = 200, n_theta = 200;

    for (int i = 0; i <= n_delta; i++)
    {
        double delta_test =
            -M_PI + i * (2 * M_PI / n_delta);

        for (int j = 0; j <= n_theta; j++)
        {
            double sin2_2theta13 =
                0.075 +
                j * (0.095 - 0.075) / n_theta;

            /* Convert to theta13 */
            double theta13_test =
                0.5 * asin(sqrt(sin2_2theta13));

            /* Projection */
            glb_projection proj =
                glbAllocProjection();

            glbDefineProjection(proj,
                                 GLB_FREE,   /* theta12 */
                                 GLB_FIXED,  /* theta13 scanned */
                                 GLB_FREE,   /* theta23 */
                                 GLB_FIXED,  /* deltaCP scanned */
                                 GLB_FREE,
                                 GLB_FREE);

            glbSetProjection(proj);

            /* Test point */
            glbCopyParams(true_values, test_values);

            glbSetOscParams(test_values,
                            theta13_test,
                            GLB_THETA_13);

            glbSetOscParams(test_values,
                            delta_test,
                            GLB_DELTA_CP);

            double chi2 =
                glbChiNP(test_values, minimum, GLB_ALL);

            double delta_chi2 =
                chi2 - chi2_min;

            fprintf(fp, "%f %f %f\n",
                    delta_test / M_PI,
                    sin2_2theta13,
                    delta_chi2);

            glbFreeProjection(proj);
        }

        fprintf(fp, "\n");
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