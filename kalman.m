clear all
graphics_toolkit("gnuplot");
figure('Position',[100, 100, 1000, 800]);

dt = 0.1;
duration = 3;
randn
w = 1 * [(dt^2/2)*randn; dt*randn];
v = 1 * [randn];
u = 2;

# Covariance
Q = [w(1)*w(1), w(1)*w(2); w(2)*w(1), w(2)*w(2)];
R = v^2;
P = Q;

# Initial
X = [0; 0];

# Model
A = [1, dt; 0, 1]; 
B = [dt^2/2; dt];
H = [1, 0];
Z = [0];

X_pred = [];
X_est = [];
Z_accum = [];

Xp = X;

for t = 0 : dt : duration

    # Save for plotting
    X_pred = [X_pred, Xp];
    X_est = [X_est, X];
    Z_accum = [Z_accum, Z]

    # we will inject this error in process
    w = 0.9 * [randn; 0];
    w = 0.01 * [(dt^2/2)*randn; dt*randn];
    err_measure = 0.8 * [randn];

    # Do KALMAN
    Xp = A*Xp + B*u + w;
    Pp = A*P*A' + Q;

    K = Pp*H' + inv(H*P*H' + R);

    Z = Xp(1,:) + err_measure;
    X = Xp + K*(Z - H*Xp);
    P = Pp - K*H*Pp;

end

tmp = X_pred;
tmp(2,:) = [];
plot(0:dt:duration, tmp, "-r");
hold on

tmp = X_est;
tmp
tmp(2,:) = [];
plot(0:dt:duration, tmp, "-g");
hold on

tmp = Z_accum;
plot(0:dt:duration, tmp, "-b");
hold on

xlabel("time (sec)");
ylabel("distance");
legend("Prediction", "Estimate", "Measurement");
pause;



