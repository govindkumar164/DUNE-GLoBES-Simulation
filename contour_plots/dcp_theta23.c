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
                     0.5878 * 0.0215,
                     0.1487 * 0.0135,
                     0.8465 * 0.0326,
                     0.0,
                     7.49e-5 * 0.0251,
                     2.534e-3 * 0.0094);

    glbSetDensityParams(input_errors, 0.05, GLB_ALL);

    /* True values */
    double delta_true = -0.79 * M_PI;
    double sin2_theta23_true = 0.561;

    /* Convert sin²θ23 → θ23 */
    double theta23_true =
        asin(sqrt(sin2_theta23_true));

    glbDefineParams(true_values,
                     0.5878,
                     0.1487,
                     theta23_true,
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
    FILE *fp = fopen("contour_theta23.dat", "w");

    /* Scan grid */
    int n_delta = 200, n_theta23 = 200;

    for (int i = 0; i <= n_delta; i++)
    {
        double delta_test =
            -M_PI + i * (2 * M_PI / n_delta);

        for (int j = 0; j <= n_theta23; j++)
        {
            /* NuFIT 3σ range */
            double sin2_theta23 =
                0.43 +
                j * (0.60 - 0.43) / n_theta23;

            /* Convert to theta23 */
            double theta23_test =
                asin(sqrt(sin2_theta23));

            /* Projection */
            glb_projection proj =
                glbAllocProjection();

            glbDefineProjection(proj,
                                 GLB_FREE,   /* theta12 */
                                 GLB_FREE,   /* theta13 */
                                 GLB_FIXED,  /* theta23 scanned */
                                 GLB_FIXED,  /* deltaCP scanned */
                                 GLB_FREE,
                                 GLB_FREE);

            glbSetProjection(proj);

            /* Test point */
            glbCopyParams(true_values, test_values);

            glbSetOscParams(test_values,
                            theta23_test,
                            GLB_THETA_23);

            glbSetOscParams(test_values,
                            delta_test,
                            GLB_DELTA_CP);

            double chi2 =
                glbChiNP(test_values, minimum, GLB_ALL);

            double delta_chi2 =
                chi2 - chi2_min;

            fprintf(fp, "%f %f %f\n",
                    delta_test / M_PI,
                    sin2_theta23,
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