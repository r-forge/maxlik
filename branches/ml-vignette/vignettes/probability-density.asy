unitsize(12mm);
defaultpen(fontsize(9));

// Axes
real x0 = -3.3;
real x1 = 3.3;
real y1 = 0.5;
path xaxis = (x0,0)--(x1,0);
path yaxis = (0,0)--(0,y1);
draw(xaxis, Arrow(TeXHead, 1));
draw(yaxis, Arrow(TeXHead, 1));
label("$x$", point(xaxis, 1), E);
label("$y$", point(yaxis, 1), SW);

// Axis labels
for(int x = (int)x0; x <= (int)x1; ++x) {
  draw((x,0)--(x,-0.1));
  label(string(x), (x,-0.1), S);
}
for(int y = 1; y < (int)y1; ++y) {
  draw((0,y)--(-0.1,y));
  label(string(y), (-0.1,y), W);
}

// normal density
real dnorm(real x) {
  1/sqrt(2*pi)*exp(-1/2*x^2);
}

path normalCurve;

for(real x = x0; x += 0.2; x < x1) {
  normalCurve = normalCurve..(x, dnorm(x));
}

// regression line
// betas: note: because of the randomness these will not be quite correct
real b0 = 0.1;
real b1 = 0.5;
real predict(real x) {
  return b0 + b1*x;
}
path regression_line = (x0, predict(x0))--(x1, predict(x1));
draw(regression_line, linewidth(1));

// Data points
srand(0);
int N=12;
real[] x;
for(int i = 0; i < N; ++i) {
  x.push(x0 + unitrand()*(x1 - x0));
}
pair[] data;
for(real xi : x) {
  real y = predict(xi) + 1.5*(unitrand() - 0.5);
  data.push((xi, y));
}
for(pair dot : data) {
  dot(dot, mediumblue);
}

// intercept
pair b = intersectionpoint(yaxis, regression_line);
draw(b--(b.x+0.1,b.y));
path intercept = (0.1,0)--(b.x+0.1,b.y);
draw(intercept, dashed);
label("$\beta_0$", point(intercept, 0.5), E);

// slope
for(int x = (int)point(regression_line, 0).x;
    x < (int)point(regression_line, 1).x;
    ++x) {
  pair prediction = intersectionpoint(regression_line, (x,0)--(x,10));
  pair cpred = (prediction.x + 1, prediction.y);
  pair prediction1 = intersectionpoint(regression_line, (x+1,0)--(x+1,10));
  draw((prediction1.x, 0)--cpred--prediction, dotted);
  path slope = prediction1--cpred;
  draw(slope, dashed);
  label("$\beta_1$", point(slope, 0.5), E);
}
