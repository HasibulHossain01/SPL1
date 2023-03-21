#include <iostream>
#include <cmath>

using namespace std;

const double g = 9.81; // acceleration due to gravity
const double rho = 1.225; // air density
const double C = 0.47; // drag coefficient
const double A = 0.05; // cross-sectional area
const double m = 0.15; // mass of projectile

double f1(double vx) {
    return -rho*C*A*vx*vx/(2.0*m);
}

double f2(double vy) {
    return -g-rho*C*A*vy*vy/(2.0*m);
}

void rk4(double &x, double &y, double &vx, double &vy, double dt) {
    double k1x = vx*dt;
    double k1y = vy*dt;
    double k1vx = f1(vx)*dt;
    double k1vy = f2(vy)*dt;

    double k2x = (vx + 0.5*k1vx)*dt;
    double k2y = (vy + 0.5*k1vy)*dt;
    double k2vx = f1(vx + 0.5*k1vx)*dt;
    double k2vy = f2(vy + 0.5*k1vy)*dt;

    double k3x = (vx + 0.5*k2vx)*dt;
    double k3y = (vy + 0.5*k2vy)*dt;
    double k3vx = f1(vx + 0.5*k2vx)*dt;
    double k3vy = f2(vy + 0.5*k2vy)*dt;

    double k4x = (vx + k3vx)*dt;
    double k4y = (vy + k3vy)*dt;
    double k4vx = f1(vx + k3vx)*dt;
    double k4vy = f2(vy + k3vy)*dt;

    x += (k1x + 2.0*k2x + 2.0*k3x + k4x)/6.0;
    y += (k1y + 2.0*k2y + 2.0*k3y + k4y)/6.0;
    vx += (k1vx + 2.0*k2vx + 2.0*k3vx + k4vx)/6.0;
    vy += (k1vy + 2.0*k2vy + 2.0*k3vy + k4vy)/6.0;
}

int main() {
    double x = 0.0, y = 0.0, vx = 30.0, vy = 30.0;
    double dt = 0.01;

    while (y >= 0.0) {
        cout << x << "," << y << endl;
        rk4(x, y, vx, vy, dt);
    }

    return 0;
}
