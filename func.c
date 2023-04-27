double t_max = 0.5;
double x_max = 0.6;
double t_step = 0.01;
double x_step = 0.01;


double func(double t, double x){
    return x*t;
}

double fi(double x){
    return x*x*x / 12;
}

double ksi(double t){
    return t*t*t / 12;
}

// return u[m, k + 1], according to sceme:
//          -if k > 0 - cross-sceme 
//          -if (k = 0) - left angle
//          -if (m = M - 1) - left angle ('cause we have no conditions for right edge)
//          -if m = 0 - from initial conditions
double *sceme_realization(double *u, int k, int m){
    int K = (int)(t_max/t_step);
    int M = (int)(x_max/x_step);
    if (m == 0){
        u[m + (k + 1) * K] = ksi( (k + 1) * t_step);
    }
    if ((k == 0) || (m == M - 1)){
        u[m + (k + 1) * K] = func(t_step * k, x_step * m) * t_step - u[k * K + m] - t_step * (u[k * K + m] - u[k * K + m - 1]) / x_step;
    }
    u[(k + 1) * K + m ] = func(t_step * k, x_step * m) * t_step + u[(k - 1) * K + m] - t_step * (u[k * K + m + 1] - u[k * K + m - 1]) / x_step;

    return u;
}

