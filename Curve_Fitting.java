package com.example;

public class Curve_Fitting {
    public void Curve_Fitting(int num, double[] xPoints, double[] yPoints, double point) {
        double sumx = 0;
        for (int i = 0; i < num; i++) {
            sumx += xPoints[i];
        }
        double sumy = 0;
        for (int i = 0; i < num; i++) {
            sumy += yPoints[i];
        }
        double sumx2 = 0;
        for (int i = 0; i < num; i++) {
            sumx2 += (xPoints[i] * xPoints[i]);
        }
        double sumxy = 0;
        for (int i = 0; i < num; i++) {
            sumxy += (xPoints[i] * yPoints[i]);
        }
        System.out.println("sum of x = " + sumx);
        System.out.println("sum of y = " + sumy);
        System.out.println("sum of x^2 = " + sumx2);
        System.out.println("sum of xy = " + sumxy);
        double a1 = sumx, b1 = num, c1 = sumy;
        double a2 = sumx2, b2 = sumx, c2 = sumxy;

        double D = a1 * b2 - a2 * b1;
        double Dx = c1 * b2 - c2 * b1;
        double Dy = a1 * c2 - a2 * c1;
        double a = 0;
        double b = 0;

        if (D != 0) {
            a = Dx / D;
            b = Dy / D;
            System.out.println("Solution:");
            System.out.println("a = " + a);
            System.out.println("b = " + b);
        } else {
            if (Dx == 0 && Dy == 0) {
                System.out.println("Infinite solutions.");
            } else {
                System.out.println("No solution.");
            }
        }
        System.out.println("y = " + a + " x + " + b);
        double result = (a * point) + b;
        System.out.println("value of y at point " + point + " = " + result);
    }
}