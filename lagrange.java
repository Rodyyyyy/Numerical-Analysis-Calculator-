package com.example;

public class lagrange {
    public void lagrange(int n, double[] x, double[] y, double xValue) {
        double result = 0;
        for (int i = 0; i < n; i++) {
            double term = y[i];
            for (int j = 0; j < n; j++) {
                if (j != i) {
                    term *= (xValue - x[j]) / (x[i] - x[j]);
                }
            }
            result += term;
        }
        System.out.println("The interpolated value at x = " + xValue + " is " + result);
    }
}