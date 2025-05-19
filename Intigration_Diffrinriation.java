package com.example;

import net.objecthunter.exp4j.Expression;
import net.objecthunter.exp4j.ExpressionBuilder;

public class Intigration_Diffrinriation {
    public void trapezoidal(String fx, double upper, double lower) {
        Expression expression1 = new ExpressionBuilder(fx)
                .variable("x")
                .build()
                .setVariable("x", upper);
        Expression expression2 = new ExpressionBuilder(fx)
                .variable("x")
                .build()
                .setVariable("x", lower);

        double result1 = expression1.evaluate();
        double result2 = expression2.evaluate();
        double result = (upper - lower) * (result1 + result2) / 2;
        System.out.println("result of the trapezoidal method is " + result);
    }

    public void composite_trapezoidal(String fx, double upper, double lower, int n) {
        double h = (upper - lower) / n;
        Expression expression1 = new ExpressionBuilder(fx)
                .variable("x")
                .build()
                .setVariable("x", upper);
        Expression expression2 = new ExpressionBuilder(fx)
                .variable("x")
                .build()
                .setVariable("x", lower);

        double result1 = expression1.evaluate();
        double result2 = expression2.evaluate();
        double sum = 0;
        for (int i = 1; i < n; i++) {
            if (i * h + lower < upper) {
                Expression expression3 = new ExpressionBuilder(fx)
                        .variable("x")
                        .build()
                        .setVariable("x", lower + i * h);
                double result3 = expression3.evaluate();
                sum += result3;
            }
        }
        double result = (h / 2) * (result1 + result2 + 2 * sum);
        System.out.println("result of the composite trapezoidal method is " + result);
    }

    public void simpson(String fx, double upper, double lower) {
        Expression expression1 = new ExpressionBuilder(fx)
                .variable("x")
                .build()
                .setVariable("x", upper);
        Expression expression2 = new ExpressionBuilder(fx)
                .variable("x")
                .build()
                .setVariable("x", lower);

        double result1 = expression1.evaluate();
        double result2 = expression2.evaluate();
        double result = (upper - lower) * (result1 + 4 * result2) / 6;
        System.out.println("result of the simpson method is " + result);
    }

    public void composite_simpson(String fx, double upper, double lower, int n) {
        double h = (upper - lower) / n;
        Expression expression1 = new ExpressionBuilder(fx)
                .variable("x")
                .build()
                .setVariable("x", upper);
        Expression expression2 = new ExpressionBuilder(fx)
                .variable("x")
                .build()
                .setVariable("x", lower);

        double result1 = expression1.evaluate();
        double result2 = expression2.evaluate();
        double sum1 = 0;
        double sum2 = 0;
        for (int i = 1; i < n; i++) {
            if (i * h + lower < upper) {
                Expression expression3 = new ExpressionBuilder(fx)
                        .variable("x")
                        .build()
                        .setVariable("x", lower + i * h);
                double result3 = expression3.evaluate();
                if (i % 2 == 0) {
                    sum1 += result3;
                } else {
                    sum2 += result3;
                }
            }
        }
        double result = (h / 3) * (result1 + 4 * sum2 + 2 * sum1 + result2);
        System.out.println("result of the composite simpson method is " + result);
    }

    public void two_points_forward(int num, double[] pointsx, double[] pointsy, double point) {
        for (int i = 0; i < num; i++) {
            if (pointsx[i] == point) {
                if (i == num - 1) {
                    System.out.println("the point is not in the range");
                }
            }
        }
        double h = pointsx[1] - pointsx[0];
        for (int i = 0; i < num; i++) {
            if (pointsx[i] == point) {
                double result = (pointsy[i + 1] - pointsy[i]) / h;
                System.out.println("the derivative at point " + point + " is " + result);
            }
        }
    }

    public void two_points_backward(int num, double[] pointsx, double[] pointsy, double point) {
        for (int i = 0; i < num; i++) {
            if (pointsx[i] == point) {
                if (i == 0) {
                    System.out.println("the point is not in the range");
                }
            }
        }
        double h = pointsx[1] - pointsx[0];
        for (int i = 0; i < num; i++) {
            if (pointsx[i] == point) {
                double result = (pointsy[i] - pointsy[i - 1]) / h;
                System.out.println("the derivative at point " + point + " is " + result);
            }
        }
    }

    public void three_points_forward(int num, double[] pointsx, double[] pointsy, double point) {
        double h = pointsx[1] - pointsx[0];
        for (int i = 0; i < num; i++) {
            if (pointsx[i] == point) {
                if (i == num - 1 || i == num - 2) {
                    System.out.println("the point is not in the range");
                }
            }
        }
        for (int i = 0; i < num; i++) {
            if (pointsx[i] == point) {
                double result = (pointsy[i] * -3 + 4 * pointsy[i + 1] - pointsy[i + 2]) / (2 * h);
                System.out.println("the derivative at point " + point + " is " + result);
            }
        }
    }

    public void three_points_backward(int num, double[] pointsx, double[] pointsy, double point) {
        for (int i = 0; i < num; i++) {
            if (pointsx[i] == point) {
                if (i == 0 || i == 1) {
                    System.out.println("the point is not in the range");
                }
            }
        }
        double h = pointsx[1] - pointsx[0];
        for (int i = 0; i < num; i++) {
            if (pointsx[i] == point) {
                double result = (3 * pointsy[i] - 4 * pointsy[i - 1] + pointsy[i - 2]) / (2 * h);
                System.out.println("the derivative at point " + point + " is " + result);
            }
        }
    }

    public void central(int num, double[] pointsx, double[] pointsy, double point) {
        for (int i = 0; i < num; i++) {
            if (pointsx[i] == point) {
                if (i == 0 || i == num - 1) {
                    System.out.println("the point is not in the range");
                    return;
                }
            }
        }
        double h = pointsx[1] - pointsx[0];
        for (int i = 0; i < num; i++) {
            if (pointsx[i] == point) {
                double result = (pointsy[i + 1] - pointsy[i - 1]) / (2 * h);
                System.out.println("the derivative at point " + point + " is " + result);
            }
        }
    }

    public void romberg(String fx, double upper, double lower, int n) {
        double h = (upper - lower) / n;
        Expression expression1 = new ExpressionBuilder(fx)
                .variable("x")
                .build()
                .setVariable("x", upper);
        Expression expression2 = new ExpressionBuilder(fx)
                .variable("x")
                .build()
                .setVariable("x", lower);

        double result1 = expression1.evaluate();
        double result2 = expression2.evaluate();
        double sum = 0;
        for (int i = 1; i <= n; i++) {
            if (i * h + lower < upper) {
                Expression expression3 = new ExpressionBuilder(fx)
                        .variable("x")
                        .build()
                        .setVariable("x", lower + i * h);
                double result3 = expression3.evaluate();
                sum += result3;
            }
        }
        double result = (h / 2) * (result1 + result2 + 2 * sum);
        System.out.println("result of the romberg method is " + result);
    }

    public void gauss(String function, double a, double b) {
        double x1 = 0;
        double w1 = 2.0;
        double xi1 = ((b - a) / 2) * x1 + (a + b) / 2;
        double f1 = new ExpressionBuilder(function).variable("x").build().setVariable("x", xi1).evaluate();
        double result1 = (b - a) / 2 * w1 * f1;

        double[] x2 = { -1 / Math.sqrt(3), 1 / Math.sqrt(3) };
        double[] w2 = { 1.0, 1.0 };
        double result2 = 0.0;
        for (int i = 0; i < 2; i++) {
            double xi = ((b - a) / 2) * x2[i] + (a + b) / 2;
            double fx = new ExpressionBuilder(function).variable("x").build().setVariable("x", xi).evaluate();
            result2 += w2[i] * fx;
        }
        result2 *= (b - a) / 2;

        double[] x3 = { -Math.sqrt(3.0 / 5), 0.0, Math.sqrt(3.0 / 5) };
        double[] w3 = { 5.0 / 9, 8.0 / 9, 5.0 / 9 };
        double result3 = 0.0;
        for (int i = 0; i < 3; i++) {
            double xi = ((b - a) / 2) * x3[i] + (a + b) / 2;
            double fx = new ExpressionBuilder(function).variable("x").build().setVariable("x", xi).evaluate();
            result3 += w3[i] * fx;
        }
        result3 *= (b - a) / 2;

        System.out.println("\nApproximated Integrals:");
        System.out.printf("Gauss 1-point: %.10f%n", result1);
        System.out.printf("Gauss 2-point: %.10f%n", result2);
        System.out.printf("Gauss 3-point: %.10f%n", result3);
    }
}