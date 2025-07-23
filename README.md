# Numerical Analysis Calculator

The Numerical Analysis Calculator is a desktop application developed in Java that implements a wide range of fundamental numerical methods. It offers a graphical user interface (GUI) where users can input values, select computation methods, and instantly view results. This project is aimed at helping students, educators, and engineering professionals understand and apply numerical techniques in a user-friendly environment.

---

## Features

### Root-Finding Algorithms
- **Bisection Method**: Locates the root of a function by repeatedly dividing the interval and selecting subintervals.
- **Newton-Raphson Method**: Iterative method using function derivatives to quickly converge to a solution.
- **Secant Method**: Uses two prior approximations to estimate the next root value, ideal when derivative is not available.

### Solving Linear Systems
- **Gaussian Elimination**: Step-by-step matrix-based method to solve systems of equations using row operations.

### Numerical Integration
- **Trapezoidal Rule**: Approximates the integral of a function using trapezoids.
- **Simpson’s Rule**: Provides a more accurate integration result by using parabolic arcs instead of straight lines.

### Differentiation
- **Finite Difference Approximations**: Derivatives are estimated using surrounding function values.

### Interpolation and Curve Fitting
- **Lagrange Interpolation**: Estimates intermediate values by constructing a polynomial that passes through a set of known data points.
- **Polynomial Regression**: Fits a polynomial equation to a given dataset using the least squares method.

### Additional UI Features
- **Dark and Light Themes**: The interface includes custom CSS themes that users can toggle for comfortable viewing.
- **Clear Input Fields and Error Messages**: Improves usability by providing feedback on invalid inputs or calculations that may not converge.

---

## Technologies Used

- **Java** – Core language used to implement all numerical logic and structure
- **JavaFX** – Provides the GUI components and scene layout
- **CSS** – Used to style both the light and dark themes

---

## Getting Started

### Prerequisites
- Java 8 or higher
- JavaFX SDK (if not included with your JDK)
- IDE such as IntelliJ IDEA, Eclipse, or NetBeans

### Setup Instructions

1. Clone the repository:
   ```bash
   git clone https://github.com/Rodyyyyy/Numerical-Analysis-Calculator-.git
   cd Numerical-Analysis-Calculator-
