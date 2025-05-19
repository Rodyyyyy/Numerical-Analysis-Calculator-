package com.example;

import net.objecthunter.exp4j.Expression;
import net.objecthunter.exp4j.ExpressionBuilder;

public class Non_Linear {
    public void Fixed_Point(String fx, double initial, int maxIterations) {
        double previous = initial;
        Expression expression = new ExpressionBuilder(fx)
                .variable("x")
                .build()
                .setVariable("x", previous);

        double result = expression.evaluate();
        System.out.println("result of iteration num 0 :" + result);

        int iteration = 1;
        while (iteration <= maxIterations) {
            Expression expression2 = new ExpressionBuilder(fx)
                    .variable("x")
                    .build()
                    .setVariable("x", result);

            double newResult = expression2.evaluate();
            System.out.println("result of iteration num " + iteration + " : " + newResult);

            if (Math.abs(newResult - result) < 1e-6) {
                System.out.println("Converged after " + iteration + " iterations");
                System.out.println("final result is " + newResult);
                return;
            }

            result = newResult;
            iteration++;
        }
        System.out.println("final result is " + result + " (max iterations reached)");
    }

    public void bisection(String fx, double a, double b, int maxIterations) {
        Expression expression = new ExpressionBuilder(fx)
                .variable("x")
                .build();

        expression.setVariable("x", a);
        double fa = expression.evaluate();
        expression.setVariable("x", b);
        double fb = expression.evaluate();

        if (fa * fb >= 0) {
            System.out.println("f(a) and f(b) must have opposite signs.");
            return;
        }

        double c = a;
        double fc;
        int iteration = 0;

        while (iteration < maxIterations) {
            c = (a + b) / 2;
            expression.setVariable("x", c);
            fc = expression.evaluate();

            System.out.println("Iteration " + iteration + ": c = " + c + ", f(c) = " + fc);

            if (Math.abs(fc) < 1e-6) {
                System.out.println("Converged after " + iteration + " iterations");
                System.out.println("Final root approximation: x = " + c);
                return;
            }

            expression.setVariable("x", a);
            fa = expression.evaluate();

            if (fa * fc < 0) {
                b = c;
            } else {
                a = c;
            }

            iteration++;
        }

        System.out.println("Final root approximation: x = " + c + " (max iterations reached)");
    }

    public void secant(String fx, double first, double second, int maxIterations) {
        int iteration = 1;
        double result = 0;

        while (iteration <= maxIterations) {
            Expression expression1 = new ExpressionBuilder(fx)
                    .variable("x")
                    .build()
                    .setVariable("x", first);

            double result1 = expression1.evaluate();
            Expression expression2 = new ExpressionBuilder(fx)
                    .variable("x")
                    .build()
                    .setVariable("x", second);

            double result2 = expression2.evaluate();
            double result3 = second - (((second - first) * result2) / (result2 - result1));

            System.out.println("result of iteration " + iteration + ":" + result3);

            if (Math.abs(result3 - second) < 1e-6) {
                System.out.println("Converged after " + iteration + " iterations");
                System.out.println("final result = " + result3);
                return;
            }

            iteration++;
            result = result3;
            first = second;
            second = result3;
        }
        System.out.println("final result = " + result + " (max iterations reached)");
    }

    public void false_position(String fx, double a, double b, int maxIterations) {
        int iteration = 1;
        double result = 0;

        while (iteration <= maxIterations) {
            Expression expression1 = new ExpressionBuilder(fx)
                    .variable("x")
                    .build()
                    .setVariable("x", b);

            double result1 = expression1.evaluate();
            Expression expression2 = new ExpressionBuilder(fx)
                    .variable("x")
                    .build()
                    .setVariable("x", a);

            double result2 = expression2.evaluate();
            double c = b - (((b - a) * result1) / (result1 - result2));
            Expression expression3 = new ExpressionBuilder(fx)
                    .variable("x")
                    .build()
                    .setVariable("x", c);

            double result3 = expression3.evaluate();
            System.out.println("result of iteration num " + iteration + ":" + result3);

            if (Math.abs(result3) < 1e-6) {
                System.out.println("Converged after " + iteration + " iterations");
                System.out.println("final result =" + result3);
                return;
            }

            result = result3;
            if (result3 > 0) {
                b = c;
            } else {
                a = c;
            }
            iteration++;
        }
        System.out.println("final result =" + result + " (max iterations reached)");
    }

    public void newton(String function, double initialGuess, double tolerance, int maxIterations) {
        Expression expr = new ExpressionBuilder(function)
                .variable("x")
                .build();

        double x = initialGuess;

        for (int i = 0; i < maxIterations; i++) {
            expr.setVariable("x", x);
            double fx = expr.evaluate();
            double dfx = derivative(expr, x);

            if (Math.abs(dfx) < 1e-10) {
                System.out.println("Derivative near zero — might diverge.");
                break;
            }

            double xNew = x - fx / dfx;

            System.out.println("result of iteration num " + i + ":" + xNew);

            if (Math.abs(xNew - x) < tolerance) {
                System.out.println("Converged after " + (i + 1) + " iterations.");
                break;
            }

            x = xNew;
        }

        System.out.println("final result is " + x);
    }

    public static double secondDerivative(Expression expr, double x) {
        double h = 1e-5;
        expr.setVariable("x", x + h);
        double fxh = expr.evaluate();
        expr.setVariable("x", x);
        double fx = expr.evaluate();
        expr.setVariable("x", x - h);
        double fxmh = expr.evaluate();
        return (fxh - 2 * fx + fxmh) / (h * h);
    }

    public static double derivative(Expression expr, double x) {
        double h = 1e-5;
        expr.setVariable("x", x + h);
        double fxh = expr.evaluate();
        expr.setVariable("x", x);
        double fx = expr.evaluate();
        return (fxh - fx) / h;
    }

    public void halley(String function, double initialGuess, double tolerance, int maxIterations) {
        Expression expr = new ExpressionBuilder(function)
                .variable("x")
                .build();

        double x = initialGuess;

        for (int i = 0; i < maxIterations; i++) {
            expr.setVariable("x", x);
            double fx = expr.evaluate();
            double fpx = derivative(expr, x);
            double fppx = secondDerivative(expr, x);

            double denominator = 2 * fpx * fpx - fx * fppx;

            if (Math.abs(denominator) < 1e-10) {
                System.out.println("Denominator near zero — Halley might diverge.");
                break;
            }

            double xNew = x - (2 * fx * fpx) / denominator;

            System.out.println("result of iteration num " + i + ": " + xNew);

            if (Math.abs(xNew - x) < tolerance) {
                System.out.println("Converged after " + (i + 1) + " iterations.");
                break;
            }

            x = xNew;
        }

        System.out.println("Final result is " + x);
    }

    public void aitkens(String fx, double initial, double tolerance, int maxIterations) {
        Expression expr = new ExpressionBuilder(fx)
                .variable("x")
                .build();

        double x0 = initial;
        expr.setVariable("x", x0);
        double x1 = expr.evaluate();
        expr.setVariable("x", x1);
        double x2 = expr.evaluate();

        System.out.println("Initial: x0 = " + x0);
        System.out.println("x1 = " + x1);
        System.out.println("x2 = " + x2);

        double result = x0;
        for (int i = 0; i < maxIterations; i++) {
            // Compute the Aitken's delta-squared approximation
            double delta1 = x1 - x0;
            double delta2 = x2 - x1;
            double denominator = delta2 - delta1;

            if (Math.abs(denominator) < 1e-10) {
                System.out.println("Denominator near zero — Aitken's method might diverge.");
                break;
            }

            double aitkenResult = x2 - (delta2 * delta2) / denominator;
            System.out.println("result of iteration num " + i + ": " + aitkenResult);

            // Check for convergence
            if (Math.abs(aitkenResult - result) < tolerance) {
                System.out.println("Converged after " + (i + 1) + " iterations");
                System.out.println("final result = " + aitkenResult);
                return;
            }

            // Update for the next iteration
            result = aitkenResult;
            x0 = x1;
            x1 = x2;
            expr.setVariable("x", x2);
            x2 = expr.evaluate();
        }

        System.out.println("final result = " + result + " (max iterations reached)");
    }
}