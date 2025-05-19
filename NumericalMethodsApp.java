package com.example;

import java.io.IOException;
import java.util.logging.FileHandler;
import java.util.logging.Logger;
import java.util.logging.SimpleFormatter;

import javafx.animation.FadeTransition;
import javafx.animation.ScaleTransition;
import javafx.application.Application;
import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.Scene;
import javafx.scene.control.Alert;
import javafx.scene.control.Button;
import javafx.scene.control.ChoiceBox;
import javafx.scene.control.Label;
import javafx.scene.control.ProgressIndicator;
import javafx.scene.control.ScrollPane;
import javafx.scene.control.Separator;
import javafx.scene.control.TextArea;
import javafx.scene.control.TextField;
import javafx.scene.control.ToggleButton;
import javafx.scene.control.Tooltip;
import javafx.scene.effect.DropShadow;
import javafx.scene.image.Image;
import javafx.scene.image.ImageView;
import javafx.scene.input.KeyCode;
import javafx.scene.layout.BorderPane;
import javafx.scene.layout.GridPane;
import javafx.scene.layout.HBox;
import javafx.scene.layout.Region;
import javafx.scene.layout.StackPane;
import javafx.scene.layout.VBox;
import javafx.scene.paint.Color;
import javafx.scene.text.Font;
import javafx.stage.Stage;
import javafx.util.Duration;

class System_Linear_Non {
    public void jacobi(int n, double[][] A, double[] b, double[] x, int iterations) {
        if (n <= 0 || A == null || b == null || x == null || A.length != n || b.length != n || x.length != n) {
            System.out.println("Error: Invalid input dimensions");
            return;
        }
        for (int i = 0; i < n; i++) {
            double sum = 0;
            for (int j = 0; j < n; j++) {
                if (j != i) sum += Math.abs(A[i][j]);
            }
            if (Math.abs(A[i][i]) <= sum || Math.abs(A[i][i]) < 1e-10) {
                System.out.println("Error: Matrix is not strictly diagonally dominant or has zero diagonal element at row " + (i + 1));
                return;
            }
        }
        double[] xNext = new double[n];
        System.out.println("Starting Jacobi iteration with initial guess: " + arrayToString(x));
        for (int iter = 0; iter < iterations; iter++) {
            for (int i = 0; i < n; i++) {
                double sum = 0;
                for (int j = 0; j < n; j++) {
                    if (j != i) {
                        sum += A[i][j] * x[j];
                    }
                }
                xNext[i] = (b[i] - sum) / A[i][i];
            }
            double maxDiff = 0;
            for (int i = 0; i < n; i++) {
                maxDiff = Math.max(maxDiff, Math.abs(xNext[i] - x[i]));
            }
            System.arraycopy(xNext, 0, x, 0, n);
            System.out.println("Iteration " + (iter + 1) + ": x = " + arrayToString(x) + ", Max difference: " + maxDiff);
            if (maxDiff < 1e-6) {
                System.out.println("Converged after " + (iter + 1) + " iterations.");
                break;
            }
        }
        System.out.println("\nFinal solution after " + iterations + " iterations: x = " + arrayToString(x));
    }

    private String arrayToString(double[] arr) {
        StringBuilder sb = new StringBuilder("[");
        for (int i = 0; i < arr.length; i++) {
            sb.append(String.format("%.6f", arr[i]));
            if (i < arr.length - 1) sb.append(", ");
        }
        sb.append("]");
        return sb.toString();
    }

    public void gauss(int n, double[][] A, double[] B, double[] x, int iterations) {
        for (int itr = 1; itr <= iterations; itr++) {
            for (int i = 0; i < n; i++) {
                double sum = B[i];
                for (int j = 0; j < n; j++) {
                    if (j != i) {
                        sum -= A[i][j] * x[j];
                    }
                }
                x[i] = sum / A[i][i];
            }

            System.out.print("Iteration " + itr + ": ");
            for (int i = 0; i < n; i++) {
                System.out.printf("x[%d] = %.6f ", i, x[i]);
            }
            System.out.println();
        }
    }

    public void newtonRaphson(String[] functions, double[] initial, double tolerance, int maxIterations) {
        int n = initial.length;
        double[] x = initial.clone();
        double[] xNext = new double[n];
        double[] F = new double[n];
        double[][] J = new double[n][n];
        double h = 1e-6;
        System.out.println("Starting Newton-Raphson with initial guess: " + arrayToString(initial));
        for (int iter = 0; iter < maxIterations; iter++) {
            for (int i = 0; i < n; i++) {
                F[i] = evaluateFunction(functions[i], x);
            }
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    double[] xTemp = x.clone();
                    xTemp[j] += h;
                    double fPlusH = evaluateFunction(functions[i], xTemp);
                    double f = evaluateFunction(functions[i], x);
                    J[i][j] = (fPlusH - f) / h;
                }
            }
            double[] deltaX = solveLinearSystem(J, negate(F));
            for (int i = 0; i < n; i++) {
                xNext[i] = x[i] + deltaX[i];
            }
            double maxDiff = 0;
            for (int i = 0; i < n; i++) {
                maxDiff = Math.max(maxDiff, Math.abs(xNext[i] - x[i]));
            }
            System.out.println("Iteration " + (iter + 1) + ":");
            System.out.print("x = [" + arrayToString(xNext) + "]");
            System.out.println(" Max difference: " + maxDiff);
            if (maxDiff < tolerance) {
                System.out.println("Converged after " + (iter + 1) + " iterations.");
                break;
            }
            x = xNext.clone();
        }
        System.out.println("\nFinal solution: x = [" + arrayToString(xNext) + "]");
    }

    private double evaluateFunction(String function, double[] x) {
        double result = 0;
        if (function.contains("cos")) {
            double x1 = x[0], x2 = x[1], x3 = x[2];
            result = 3 * x1 - Math.cos(x2 * x3) - 0.5;
        } else if (function.contains("625")) {
            double x1 = x[0], x2 = x[1];
            result = 4 * x1 * x1 - 625 * x2 * x2 + 2 * x2 - 1;
        } else if (function.contains("exp")) {
            double x1 = x[0], x2 = x[1], x3 = x[2];
            result = Math.exp(-x1 * x2) + 20 * x3 + (10 * Math.PI - 3) / 3;
        }
        return result;
    }

    private double[] negate(double[] vector) {
        double[] result = new double[vector.length];
        for (int i = 0; i < vector.length; i++) {
            result[i] = -vector[i];
        }
        return result;
    }

    private double[] solveLinearSystem(double[][] A, double[] b) {
        int n = b.length;
        double[][] augmented = new double[n][n + 1];
        double[] deltaX = new double[n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                augmented[i][j] = A[i][j];
            }
            augmented[i][n] = b[i];
        }
        for (int i = 0; i < n; i++) {
            double pivot = augmented[i][i];
            if (Math.abs(pivot) < 1e-10) {
                throw new RuntimeException("Matrix is singular or nearly singular");
            }
            for (int j = 0; j <= n; j++) {
                augmented[i][j] /= pivot;
            }
            for (int k = 0; k < n; k++) {
                if (k != i) {
                    double factor = augmented[k][i];
                    for (int j = 0; j <= n; j++) {
                        augmented[k][j] -= factor * augmented[i][j];
                    }
                }
            }
        }
        for (int i = 0; i < n; i++) {
            deltaX[i] = augmented[i][n];
        }
        return deltaX;
    }
}

class Non_Linear {
    public void Fixed_Point(String function, double initial, int iterations) {}
    public void bisection(String function, double initial, double second, int iterations) {}
    public void newton(String function, double initial, double tolerance, int iterations) {}
    public void secant(String function, double initial, double second, int iterations) {}
    public void false_position(String function, double initial, double second, int iterations) {}
    public void halley(String function, double initial, double tolerance, int iterations) {}
    public void aitkens(String function, double initial, double tolerance, int iterations) {}
}

class Intigration_Diffrinriation {
    public void trapezoidal(String function, double upper, double lower) {
        if (upper <= 0 || lower <= 0) {
            System.out.println("Error: ln is undefined for non-positive values. Ensure upper and lower limits are greater than 0.");
            return;
        }
        java.util.function.DoubleFunction<Double> f = x -> {
            if (function.equals("ln(x)")) {
                return Math.log(x);
            } else {
                System.out.println("Error: Unsupported function. Only 'ln(x)' is supported for now.");
                return 0.0;
            }
        };
        double a = lower;
        double b = upper;
        double fa = f.apply(a);
        double fb = f.apply(b);
        double area = (b - a) / 2.0 * (fa + fb);
        System.out.println("Trapezoidal approximation of ∫ from " + a + " to " + b + " of " + function + " = " + area);
    }

    public void composite_trapezoidal(String function, double upper, double lower, int intervals) {
        if (intervals <= 0) {
            System.out.println("Error: Number of intervals must be positive.");
            return;
        }
        double a = lower;
        double b = upper;
        double h = (b - a) / intervals;
        double sum = 0.0;

        java.util.function.DoubleFunction<Double> f = x -> evaluatePolynomial(function, x);

        double fa = f.apply(a);
        double fb = f.apply(b);
        for (int i = 1; i < intervals; i++) {
            double x = a + i * h;
            sum += f.apply(x);
        }

        double area = (h / 2.0) * (fa + 2 * sum + fb);
        System.out.println("Composite Trapezoidal approximation of ∫ from " + a + " to " + b + " of " + function + " with " + intervals + " intervals = " + area);
    }

    private double evaluatePolynomial(String function, double x) {
        String[] terms = function.replaceAll("\\s+", "").split("(?=[-+])");
        double result = 0.0;

        for (String term : terms) {
            double coefficient = 1.0;
            double exponent = 0.0;
            boolean isNegative = term.startsWith("-");
            if (isNegative || term.startsWith("+")) {
                term = term.substring(1);
            }

            if (term.isEmpty()) {
                continue;
            }

            if (term.equals("x")) {
                coefficient = isNegative ? -1.0 : 1.0;
                exponent = 1.0;
            } else if (term.contains("x^")) {
                String[] parts = term.split("x\\^");
                if (parts.length != 2) {
                    throw new IllegalArgumentException("Invalid term format: " + term);
                }
                coefficient = parts[0].isEmpty() ? 1.0 : Double.parseDouble(parts[0]);
                if (isNegative) coefficient = -coefficient;
                exponent = Double.parseDouble(parts[1]);
            } else if (term.contains("x")) {
                String coeffStr = term.replace("x", "");
                coefficient = coeffStr.isEmpty() ? 1.0 : Double.parseDouble(coeffStr);
                if (isNegative) coefficient = -coefficient;
                exponent = 1.0;
            } else {
                coefficient = Double.parseDouble(term);
                if (isNegative) coefficient = -coefficient;
                exponent = 0.0;
            }

            result += coefficient * Math.pow(x, exponent);
        }
        return result;
    }

    public void simpson(String function, double upper, double lower) {}
    public void composite_simpson(String function, double upper, double lower, int intervals) {}
    public void romberg(String function, double upper, double lower, int intervals) {}
    public void gauss(String function, double upper, double lower) {}
    public void two_points_forward(int numPoints, double[] xPoints, double[] yPoints, double point) {}
    public void two_points_backward(int numPoints, double[] xPoints, double[] yPoints, double point) {}
    public void three_points_forward(int numPoints, double[] xPoints, double[] yPoints, double point) {}
    public void three_points_backward(int numPoints, double[] xPoints, double[] yPoints, double point) {}
    public void central(int numPoints, double[] xPoints, double[] yPoints, double point) {}
}

class Curve_Fitting {
    public void Curve_Fitting(int numPoints, double[] xPoints, double[] yPoints, double xValue) {}
}

class lagrange {
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

public class NumericalMethodsApp extends Application {
    private BorderPane root;
    private VBox sidebar;
    private StackPane contentArea;
    private boolean isDarkMode = false;
    private static final Logger LOGGER = Logger.getLogger(NumericalMethodsApp.class.getName());
    private ProgressIndicator progressIndicator;

    @Override
    public void start(Stage primaryStage) {
        setupLogging();
        root = new BorderPane();
        progressIndicator = new ProgressIndicator();
        progressIndicator.setVisible(false);

        setupSidebar();
        setupContentArea();
        setupThemeToggle();

        Scene scene = new Scene(root, 1200, 800);
        applyStylesheet(scene, isDarkMode ? "/dark-theme.css" : "/light-theme.css");

        primaryStage.setTitle("Numerical Analysis Calculator");
        primaryStage.setScene(scene);
        primaryStage.show();

        scene.setOnKeyPressed(e -> {
            if (e.getCode() == KeyCode.T && e.isControlDown()) {
                toggleTheme();
            }
        });
    }

    private void setupLogging() {
        try {
            FileHandler fileHandler = new FileHandler("numerical_methods.log", true);
            fileHandler.setFormatter(new SimpleFormatter());
            LOGGER.addHandler(fileHandler);
        } catch (IOException e) {
            LOGGER.severe("Failed to setup logging: " + e.getMessage());
        }
    }

    private void applyStylesheet(Scene scene, String stylesheetPath) {
        try {
            String cssPath = getClass().getResource(stylesheetPath) != null
                ? getClass().getResource(stylesheetPath).toExternalForm()
                : null;
            if (cssPath != null) {
                scene.getStylesheets().clear();
                scene.getStylesheets().add(cssPath);
                LOGGER.info("Applied stylesheet: " + stylesheetPath);
            } else {
                LOGGER.warning("Stylesheet not found: " + stylesheetPath);
                showErrorDialog("Stylesheet Error", "Could not load stylesheet: " + stylesheetPath + ". Using default styles.");
                scene.getRoot().setStyle(
                    "-fx-font-family: 'Segoe UI', Arial, sans-serif;" +
                    "-fx-background-color: #ECF0F1;"
                );
            }
        } catch (Exception e) {
            LOGGER.severe("Error loading stylesheet: " + e.getMessage());
            showErrorDialog("Stylesheet Error", "Failed to apply stylesheet: " + e.getMessage());
        }
    }

    private void toggleTheme() {
        isDarkMode = !isDarkMode;
        Scene scene = root.getScene();
        applyStylesheet(scene, isDarkMode ? "/dark-theme.css" : "/light-theme.css");
    }

    private void showErrorDialog(String title, String message) {
        Alert alert = new Alert(Alert.AlertType.ERROR);
        alert.setTitle(title);
        alert.setHeaderText(null);
        alert.setContentText(message);
        alert.showAndWait();
    }

    private void setupSidebar() {
        sidebar = new VBox(15);
        sidebar.setPadding(new Insets(20));
        sidebar.setStyle("-fx-background-color: linear-gradient(to bottom, #2C3E50, #1A252F);");

        HBox titleBox = new HBox(10);
        ImageView icon = new ImageView(new Image("file:resources/calculator.png", 24, 24, true, true));
        Label title = new Label("Numerical Analysis Calculator");
        title.setFont(Font.font("Segoe UI", 24));
        title.setTextFill(Color.WHITE);
        title.getStyleClass().add("project-title");
        titleBox.getChildren().addAll(icon, title);
        titleBox.setAlignment(Pos.CENTER_LEFT);

        Separator separator = new Separator();
        separator.setStyle("-fx-background-color: #ECF0F1; -fx-opacity: 0.3;");

        Label operationsHeader = new Label("Operations");
        operationsHeader.setFont(Font.font("Segoe UI", 16));
        operationsHeader.setTextFill(Color.WHITE);
        operationsHeader.getStyleClass().add("operations-header");

        VBox operationsBox = new VBox(10);
        operationsBox.setPadding(new Insets(10, 0, 0, 0));

        String[] operations = {"Non Linear", "Curve Fitting", "System", "Lagrange", "Integration", "Differentiation"};
        for (String op : operations) {
            Button btn = createSidebarButton(op);
            btn.setOnAction(e -> {
                showOperationPanel(op);
                animateButtonClick(btn);
            });
            operationsBox.getChildren().add(btn);
        }

        sidebar.getChildren().addAll(titleBox, separator, operationsHeader, operationsBox);
        root.setLeft(sidebar);
    }

    private Button createSidebarButton(String text) {
        Button btn = new Button(text);
        btn.setPrefWidth(200);
        btn.getStyleClass().add("sidebar-button");
        btn.setTooltip(new Tooltip("Select " + text + " operations"));

        ScaleTransition scaleIn = new ScaleTransition(Duration.millis(200), btn);
        scaleIn.setToX(1.05);
        scaleIn.setToY(1.05);

        ScaleTransition scaleOut = new ScaleTransition(Duration.millis(200), btn);
        scaleOut.setToX(1.0);
        scaleOut.setToY(1.0);

        btn.setOnMouseEntered(e -> {
            scaleIn.playFromStart();
            btn.setEffect(new DropShadow(10, Color.rgb(0, 150, 255, 0.5)));
        });
        btn.setOnMouseExited(e -> {
            scaleOut.playFromStart();
            btn.setEffect(null);
        });

        return btn;
    }

    private void setupContentArea() {
        contentArea = new StackPane();
        contentArea.setStyle("-fx-background-color: linear-gradient(to bottom right, #ECF0F1, #DDE4E6);");

        Label welcome = new Label("Select an operation from the sidebar");
        welcome.setFont(Font.font("Segoe UI", 18));
        welcome.getStyleClass().add("welcome-label");

        ScrollPane scrollPane = new ScrollPane();
        scrollPane.setContent(welcome);
        scrollPane.setFitToWidth(true);
        scrollPane.setStyle("-fx-background: transparent; -fx-background-color: transparent;");

        contentArea.getChildren().addAll(scrollPane, progressIndicator);
        root.setCenter(contentArea);
    }

    private void setupThemeToggle() {
        ToggleButton themeToggle = new ToggleButton("Dark Mode");
        themeToggle.getStyleClass().add("theme-toggle");
        themeToggle.setOnAction(e -> toggleTheme());

        HBox topBar = new HBox(10, themeToggle);
        topBar.setPadding(new Insets(10));
        topBar.setAlignment(Pos.CENTER_RIGHT);
        root.setTop(topBar);
    }

    private void animateButtonClick(Button btn) {
        ScaleTransition scale = new ScaleTransition(Duration.millis(100), btn);
        scale.setToX(0.95);
        scale.setToY(0.95);
        scale.setAutoReverse(true);
        scale.setCycleCount(2);
        scale.play();
    }

    private void showOperationPanel(String operation) {
        contentArea.getChildren().clear();
        VBox panel = new VBox(15);
        panel.setPadding(new Insets(20));
        panel.setAlignment(Pos.TOP_LEFT);
        panel.getStyleClass().add("content-panel");

        Label header = new Label(operation);
        header.setFont(Font.font("Segoe UI", 28));
        header.getStyleClass().add("panel-header");

        ScrollPane scrollPane = new ScrollPane(panel);
        scrollPane.setFitToWidth(true);
        scrollPane.setStyle("-fx-background: transparent; -fx-background-color: transparent;");

        try {
            switch (operation.toLowerCase()) {
                case "non linear":
                    String[] nonLinearMethods = {"Fixed Point", "Bisection", "Newton", "Secant", "False Position", "Halley", "Aitkens"};
                    addMethodButtons(panel, nonLinearMethods, operation);
                    break;
                case "curve fitting":
                    addCurveFittingPanel(panel);
                    break;
                case "system":
                    addSystemPanel(panel);
                    break;
                case "lagrange":
                    addLagrangePanel(panel);
                    break;
                case "integration":
                    String[] integrationMethods = {"Trapezoidal", "Composite_Trapezoidal", "Simpson", "Composite_Simpson", "Romberg", "Gauss"};
                    addMethodButtons(panel, integrationMethods, operation);
                    break;
                case "differentiation":
                    addDifferentiationPanel(panel);
                    break;
                default:
                    throw new IllegalArgumentException("Unknown operation: " + operation);
            }
        } catch (Exception e) {
            LOGGER.severe("Error displaying operation panel: " + e.getMessage());
            showErrorDialog("Operation Error", "Failed to load " + operation + " panel: " + e.getMessage());
        }

        contentArea.getChildren().addAll(scrollPane, progressIndicator);
        animatePanelTransition(scrollPane);
    }

    private void addMethodButtons(VBox panel, String[] methods, String operation) {
        GridPane buttonGrid = new GridPane();
        buttonGrid.setHgap(10);
        buttonGrid.setVgap(10);

        for (int i = 0; i < methods.length; i++) {
            Button btn = createMethodButton(methods[i]);
            final String methodName = methods[i];
            final Button btnCopy = btn;
            btn.setOnAction(e -> {
                showMethodInputPanel(operation, methodName);
                animateButtonClick(btnCopy);
            });
            buttonGrid.add(btn, i % 3, i / 3);
        }

        panel.getChildren().add(buttonGrid);
    }

    private Button createMethodButton(String text) {
        Button btn = new Button(text);
        btn.setPrefWidth(220);
        btn.getStyleClass().add("method-button");
        btn.setTooltip(new Tooltip("Calculate using " + text + " method"));

        ScaleTransition scaleIn = new ScaleTransition(Duration.millis(200), btn);
        scaleIn.setToX(1.05);
        scaleIn.setToY(1.05);

        ScaleTransition scaleOut = new ScaleTransition(Duration.millis(200), btn);
        scaleOut.setToX(1.0);
        scaleOut.setToY(1.0);

        btn.setOnMouseEntered(e -> {
            scaleIn.playFromStart();
            btn.setEffect(new DropShadow(10, Color.rgb(0, 150, 255, 0.5)));
        });
        btn.setOnMouseExited(e -> {
            scaleOut.playFromStart();
            btn.setEffect(null);
        });

        return btn;
    }

    private void animatePanelTransition(Region panel) {
        FadeTransition fade = new FadeTransition(Duration.millis(300), panel);
        fade.setFromValue(0);
        fade.setToValue(1);
        fade.play();
    }

    private void showMethodInputPanel(String operation, String method) {
        contentArea.getChildren().clear();
        VBox panel = new VBox(15);
        panel.setPadding(new Insets(20));
        panel.setAlignment(Pos.TOP_LEFT);
        panel.getStyleClass().add("content-panel");

        Label header = new Label(method + " Method");
        header.setFont(Font.font("Segoe UI", 24));
        header.getStyleClass().add("panel-header");

        panel.getChildren().add(header);

        ScrollPane scrollPane = new ScrollPane(panel);
        scrollPane.setFitToWidth(true);
        scrollPane.setStyle("-fx-background: transparent; -fx-background-color: transparent;");

        try {
            switch (operation.toLowerCase()) {
                case "non linear":
                    addNonLinearInputPanel(panel, method);
                    break;
                case "integration":
                    addIntegrationInputPanel(panel, method);
                    break;
                default:
                    throw new IllegalArgumentException("Unsupported operation: " + operation);
            }
        } catch (Exception e) {
            LOGGER.severe("Error displaying method panel: " + e.getMessage());
            showErrorDialog("Panel Error", "Failed to load " + method + " panel: " + e.getMessage());
        }

        contentArea.getChildren().addAll(scrollPane, progressIndicator);
        animatePanelTransition(scrollPane);
    }

    private void addNonLinearInputPanel(VBox panel, String method) {
        Non_Linear nonLinear = new Non_Linear();
        VBox inputBox = createInputBox();

        TextField functionField = createTextField("Enter function (e.g., x^2 - 4)");
        TextField initialField = createTextField("Enter initial guess");
        TextField secondField = createTextField("Enter second point");
        TextField toleranceField = createTextField("Enter tolerance (e.g., 1e-6)");
        TextField iterationsField = createTextField("Enter max iterations");
        TextArea outputArea = createOutputArea();

        addLabeledField(inputBox, "f(x):", functionField, "Mathematical function in terms of x");

        String initialLabel = method.equalsIgnoreCase("bisection") || method.equalsIgnoreCase("false position") || method.equalsIgnoreCase("secant")
            ? "First Point (a):" : "Initial Guess:";
        addLabeledField(inputBox, initialLabel, initialField, "Starting point for iteration");

        if (method.equalsIgnoreCase("bisection") || method.equalsIgnoreCase("false position") || method.equalsIgnoreCase("secant")) {
            String secondLabel = method.equalsIgnoreCase("secant") ? "Second Point:" : "Second Point (b):";
            addLabeledField(inputBox, secondLabel, secondField, "Second point for interval");
        }

        if (method.equalsIgnoreCase("newton") || method.equalsIgnoreCase("halley") || method.equalsIgnoreCase("aitkens")) {
            addLabeledField(inputBox, "Tolerance:", toleranceField, "Convergence tolerance");
        }

        addLabeledField(inputBox, "Max Iterations:", iterationsField, "Maximum number of iterations");

        Button calculateBtn = createMethodButton("Calculate");
        calculateBtn.setOnAction(e -> {
            progressIndicator.setVisible(true);
            try {
                validateNonLinearInputs(functionField, initialField, secondField, toleranceField, iterationsField, method);
                outputArea.clear();
                String function = functionField.getText();
                double initial = initialField.getText().isEmpty() ? 0 : Double.parseDouble(initialField.getText());
                double second = secondField.getText().isEmpty() ? 0 : Double.parseDouble(secondField.getText());
                double tolerance = toleranceField.getText().isEmpty() ? 1e-6 : Double.parseDouble(toleranceField.getText());
                int iterations = iterationsField.getText().isEmpty() ? 100 : Integer.parseInt(iterationsField.getText());

                TextAreaOutputStream outStream = new TextAreaOutputStream(outputArea);
                System.setOut(new java.io.PrintStream(outStream));

                switch (method.toLowerCase()) {
                    case "fixed point":
                        nonLinear.Fixed_Point(function, initial, iterations);
                        break;
                    case "bisection":
                        nonLinear.bisection(function, initial, second, iterations);
                        break;
                    case "newton":
                        nonLinear.newton(function, initial, tolerance, iterations);
                        break;
                    case "secant":
                        nonLinear.secant(function, initial, second, iterations);
                        break;
                    case "false position":
                        nonLinear.false_position(function, initial, second, iterations);
                        break;
                    case "halley":
                        nonLinear.halley(function, initial, tolerance, iterations);
                        break;
                    case "aitkens":
                        nonLinear.aitkens(function, initial, tolerance, iterations);
                        break;
                    default:
                        throw new IllegalArgumentException("Unknown method: " + method);
                }
            } catch (Exception ex) {
                LOGGER.warning("Calculation error: " + ex.getMessage());
                showErrorDialog("Calculation Error", getUserFriendlyErrorMessage(ex));
                outputArea.setText("Error: " + getUserFriendlyErrorMessage(ex));
            } finally {
                progressIndicator.setVisible(false);
            }
        });

        addLabeledField(inputBox, "Output:", outputArea, "Calculation results");
        inputBox.getChildren().add(calculateBtn);
        panel.getChildren().add(inputBox);
    }

    private void addSystemPanel(VBox panel) {
        System_Linear_Non system = new System_Linear_Non();
        VBox inputBox = createInputBox();

        ChoiceBox<String> typeChoice = createChoiceBox(new String[]{"Linear", "Non-Linear"}, "Select system type");
        TextField equationsField = createTextField("Enter number of equations");
        TextArea coefficientsInput = createTextArea("Enter coefficients of A (one row per line, space-separated)\nThen enter b values (space-separated)");
        TextField initialGuessesField = createTextField("Enter initial guesses (space-separated)");
        TextField toleranceField = createTextField("Enter tolerance (e.g., 1e-6)");
        TextField iterationsField = createTextField("Enter iterations");
        ChoiceBox<String> methodChoice = createChoiceBox(new String[]{}, "Select method");
        TextArea outputArea = createOutputArea();

        typeChoice.setOnAction(e -> {
            methodChoice.getItems().clear();
            if (typeChoice.getValue() != null && typeChoice.getValue().equals("Linear")) {
                methodChoice.getItems().addAll("Jacobi", "Gauss");
                coefficientsInput.setPromptText("Enter coefficients of A (one row per line, space-separated)\nThen enter b values (space-separated)");
            } else {
                methodChoice.getItems().add("Newton");
                coefficientsInput.setPromptText("Enter functions F1, F2, F3 (one per line)");
            }
        });

        addLabeledField(inputBox, "System Type:", typeChoice, "Linear or Non-Linear system");
        addLabeledField(inputBox, "Number of Equations:", equationsField, "Number of equations in system");
        addLabeledField(inputBox, "Input (A and b or Functions):", coefficientsInput, "Enter coefficients for linear systems or functions for non-linear");
        addLabeledField(inputBox, "Initial Guesses:", initialGuessesField, "Space-separated initial guesses");
        addLabeledField(inputBox, "Tolerance:", toleranceField, "Convergence tolerance");
        addLabeledField(inputBox, "Iterations:", iterationsField, "Number of iterations");
        addLabeledField(inputBox, "Method:", methodChoice, "Calculation method");

        Button calculateBtn = createMethodButton("Calculate");
        calculateBtn.setOnAction(e -> {
            progressIndicator.setVisible(true);
            try {
                validateSystemInputs(typeChoice, equationsField, coefficientsInput, initialGuessesField, toleranceField, iterationsField, methodChoice);
                outputArea.clear();
                String type = typeChoice.getValue();
                int n = Integer.parseInt(equationsField.getText());
                String[] guesses = initialGuessesField.getText().trim().split("\\s+");
                if (guesses.length != n) throw new IllegalArgumentException("Number of initial guesses must match number of equations");
                double[] initial = new double[n];
                for (int i = 0; i < n; i++) {
                    initial[i] = Double.parseDouble(guesses[i]);
                }
                double tolerance = toleranceField.getText().isEmpty() ? 1e-6 : Double.parseDouble(toleranceField.getText());
                int iterations = iterationsField.getText().isEmpty() ? 100 : Integer.parseInt(iterationsField.getText());

                TextAreaOutputStream outStream = new TextAreaOutputStream(outputArea);
                System.setOut(new java.io.PrintStream(outStream));

                if (type.equals("Linear")) {
                    String[] lines = coefficientsInput.getText().split("\n");
                    if (lines.length < n + 1) throw new IllegalArgumentException("Insufficient input: provide " + n + " rows for A and 1 row for b");
                    double[][] A = new double[n][n];
                    double[] b = new double[n];
                    for (int i = 0; i < n; i++) {
                        String[] parts = lines[i].trim().split("\\s+");
                        if (parts.length != n) throw new IllegalArgumentException("Invalid row format for A at line " + (i + 1));
                        for (int j = 0; j < n; j++) {
                            A[i][j] = Double.parseDouble(parts[j]);
                        }
                    }
                    String[] bParts = lines[n].trim().split("\\s+");
                    if (bParts.length != n) throw new IllegalArgumentException("Invalid format for b");
                    for (int i = 0; i < n; i++) {
                        b[i] = Double.parseDouble(bParts[i]);
                    }

                    if (methodChoice.getValue().equals("Jacobi")) {
                        system.jacobi(n, A, b, initial, iterations);
                    } else {
                        system.gauss(n, A, b, initial, iterations);
                    }
                } else {
                    String[] functions = coefficientsInput.getText().split("\n");
                    if (functions.length != n) throw new IllegalArgumentException("Number of functions must match number of equations");
                    if (methodChoice.getValue().equals("Newton")) {
                        system.newtonRaphson(functions, initial, tolerance, iterations);
                    } else {
                        outputArea.setText("Method not supported for non-linear systems.");
                    }
                }
            } catch (Exception ex) {
                LOGGER.warning("Calculation error: " + ex.getMessage());
                showErrorDialog("Calculation Error", getUserFriendlyErrorMessage(ex));
                outputArea.setText("Error: " + getUserFriendlyErrorMessage(ex));
            } finally {
                progressIndicator.setVisible(false);
            }
        });

        addLabeledField(inputBox, "Output:", outputArea, "System solution results");
        inputBox.getChildren().add(calculateBtn);
        panel.getChildren().add(inputBox);
    }

    private void addIntegrationInputPanel(VBox panel, String method) {
        Intigration_Diffrinriation integration = new Intigration_Diffrinriation();
        VBox inputBox = createInputBox();

        TextField functionField = createTextField("Enter function (e.g., ln(x))");
        TextField upperField = createTextField("Enter upper limit");
        TextField lowerField = createTextField("Enter lower limit");
        TextField intervalsField = createTextField("Enter number of subintervals");
        TextArea outputArea = createOutputArea();

        addLabeledField(inputBox, "f(x):", functionField, "Function to integrate (e.g., ln(x))");
        addLabeledField(inputBox, "Upper Limit:", upperField, "Upper bound of integration");
        addLabeledField(inputBox, "Lower Limit:", lowerField, "Lower bound of integration");

        if (method.equalsIgnoreCase("composite_trapezoidal") || method.equalsIgnoreCase("composite_simpson") || method.equalsIgnoreCase("romberg")) {
            addLabeledField(inputBox, "Number of Subintervals:", intervalsField, "Number of integration intervals");
        }

        Button calculateBtn = createMethodButton("Calculate");
        calculateBtn.setOnAction(e -> {
            progressIndicator.setVisible(true);
            try {
                validateIntegrationInputs(functionField, upperField, lowerField, intervalsField, method);
                outputArea.clear();
                String function = functionField.getText();
                double upper = Double.parseDouble(upperField.getText());
                double lower = Double.parseDouble(lowerField.getText());
                int intervals = intervalsField.getText().isEmpty() ? 0 : Integer.parseInt(intervalsField.getText());

                TextAreaOutputStream outStream = new TextAreaOutputStream(outputArea);
                System.setOut(new java.io.PrintStream(outStream));

                switch (method.toLowerCase()) {
                    case "trapezoidal":
                        integration.trapezoidal(function, upper, lower);
                        break;
                    case "composite_trapezoidal":
                        integration.composite_trapezoidal(function, upper, lower, intervals);
                        break;
                    case "simpson":
                        integration.simpson(function, upper, lower);
                        break;
                    case "composite_simpson":
                        integration.composite_simpson(function, upper, lower, intervals);
                        break;
                    case "romberg":
                        integration.romberg(function, upper, lower, intervals);
                        break;
                    case "gauss":
                        integration.gauss(function, upper, lower);
                        break;
                    default:
                        throw new IllegalArgumentException("Unknown method: " + method);
                }
            } catch (Exception ex) {
                LOGGER.warning("Calculation error: " + ex.getMessage());
                showErrorDialog("Calculation Error", getUserFriendlyErrorMessage(ex));
                outputArea.setText("Error: " + getUserFriendlyErrorMessage(ex));
            } finally {
                progressIndicator.setVisible(false);
            }
        });

        addLabeledField(inputBox, "Output:", outputArea, "Integration results");
        inputBox.getChildren().add(calculateBtn);
        panel.getChildren().add(inputBox);
    }

    private void addCurveFittingPanel(VBox panel) {
        Curve_Fitting curveFitting = new Curve_Fitting();
        VBox inputBox = createInputBox();

        TextField pointsField = createTextField("Enter number of points");
        TextArea pointsInput = createTextArea("Enter points (x y per line)");
        TextField xValueField = createTextField("Enter x value to predict");
        TextArea outputArea = createOutputArea();

        addLabeledField(inputBox, "Number of Points:", pointsField, "Total number of data points");
        addLabeledField(inputBox, "Enter Points (x y per line):", pointsInput, "Enter x and y coordinates");
        addLabeledField(inputBox, "X Value to Predict:", xValueField, "X value for prediction");

        Button calculateBtn = createMethodButton("Calculate");
        calculateBtn.setOnAction(e -> {
            progressIndicator.setVisible(true);
            try {
                validateCurveFittingInputs(pointsField, pointsInput, xValueField);
                outputArea.clear();
                int numPoints = Integer.parseInt(pointsField.getText());
                String[] lines = pointsInput.getText().split("\n");
                double[] xPoints = new double[numPoints];
                double[] yPoints = new double[numPoints];
                for (int i = 0; i < numPoints; i++) {
                    String[] parts = lines[i].trim().split("\\s+");
                    if (parts.length != 2) throw new IllegalArgumentException("Invalid point format at line " + (i + 1));
                    xPoints[i] = Double.parseDouble(parts[0]);
                    yPoints[i] = Double.parseDouble(parts[1]);
                }
                double xValue = Double.parseDouble(xValueField.getText());

                TextAreaOutputStream outStream = new TextAreaOutputStream(outputArea);
                System.setOut(new java.io.PrintStream(outStream));

                curveFitting.Curve_Fitting(numPoints, xPoints, yPoints, xValue);
            } catch (Exception ex) {
                LOGGER.warning("Calculation error: " + ex.getMessage());
                showErrorDialog("Calculation Error", getUserFriendlyErrorMessage(ex));
                outputArea.setText("Error: " + getUserFriendlyErrorMessage(ex));
            } finally {
                progressIndicator.setVisible(false);
            }
        });

        addLabeledField(inputBox, "Output:", outputArea, "Curve fitting results");
        inputBox.getChildren().add(calculateBtn);
        panel.getChildren().add(inputBox);
    }

    private void addLagrangePanel(VBox panel) {
        lagrange lagrange = new lagrange();
        VBox inputBox = createInputBox();

        TextField pointsField = createTextField("Enter number of points");
        TextArea pointsInput = createTextArea("Enter points (x y per line)");
        TextField xValueField = createTextField("Enter x value to interpolate");
        TextArea outputArea = createOutputArea();

        addLabeledField(inputBox, "Number of Points:", pointsField, "Total number of data points");
        addLabeledField(inputBox, "Enter Points (x y per line):", pointsInput, "Enter x and y coordinates");
        addLabeledField(inputBox, "X Value to Interpolate:", xValueField, "X value for interpolation");

        Button calculateBtn = createMethodButton("Calculate");
        calculateBtn.setOnAction(e -> {
            progressIndicator.setVisible(true);
            try {
                validateLagrangeInputs(pointsField, pointsInput, xValueField);
                outputArea.clear();
                int n = Integer.parseInt(pointsField.getText());
                String[] lines = pointsInput.getText().split("\n");
                double[] x = new double[n];
                double[] y = new double[n];
                for (int i = 0; i < n; i++) {
                    String[] parts = lines[i].trim().split("\\s+");
                    if (parts.length != 2) throw new IllegalArgumentException("Invalid point format at line " + (i + 1));
                    x[i] = Double.parseDouble(parts[0]);
                    y[i] = Double.parseDouble(parts[1]);
                }
                double xValue = Double.parseDouble(xValueField.getText());

                TextAreaOutputStream outStream = new TextAreaOutputStream(outputArea);
                System.setOut(new java.io.PrintStream(outStream));

                lagrange.lagrange(n, x, y, xValue);
            } catch (Exception ex) {
                LOGGER.warning("Calculation error: " + ex.getMessage());
                showErrorDialog("Calculation Error", getUserFriendlyErrorMessage(ex));
                outputArea.setText("Error: " + getUserFriendlyErrorMessage(ex));
            } finally {
                progressIndicator.setVisible(false);
            }
        });

        addLabeledField(inputBox, "Output:", outputArea, "Interpolation results");
        inputBox.getChildren().add(calculateBtn);
        panel.getChildren().add(inputBox);
    }

    private void addDifferentiationPanel(VBox panel) {
        Intigration_Diffrinriation differentiation = new Intigration_Diffrinriation();
        VBox inputBox = createInputBox();

        ChoiceBox<String> pointsChoice = createChoiceBox(new String[]{"2 Points", "3 Points"}, "Select number of points");
        ChoiceBox<String> methodChoice = createChoiceBox(new String[]{}, "Select differentiation method");
        TextField pointsField = createTextField("Enter number of points");
        TextArea pointsInput = createTextArea("Enter points (x y per line)");
        TextField xValueField = createTextField("Enter point to differentiate");
        TextArea outputArea = createOutputArea();

        pointsChoice.setOnAction(e -> {
            methodChoice.getItems().clear();
            if (pointsChoice.getValue() != null && pointsChoice.getValue().equals("2 Points")) {
                methodChoice.getItems().addAll("Forward", "Backward");
            } else if (pointsChoice.getValue() != null) {
                methodChoice.getItems().addAll("Forward", "Backward", "Central");
            }
        });

        addLabeledField(inputBox, "Number of Points:", pointsChoice, "2 or 3 point method");
        addLabeledField(inputBox, "Method:", methodChoice, "Differentiation method");
        addLabeledField(inputBox, "Number of Points:", pointsField, "Total number of data points");
        addLabeledField(inputBox, "Enter Points (x y per line):", pointsInput, "Enter x and y coordinates");
        addLabeledField(inputBox, "Point to Differentiate:", xValueField, "Point for differentiation");

        Button calculateBtn = createMethodButton("Calculate");
        calculateBtn.setOnAction(e -> {
            progressIndicator.setVisible(true);
            try {
                validateDifferentiationInputs(pointsChoice, methodChoice, pointsField, pointsInput, xValueField);
                outputArea.clear();
                int numPoints = Integer.parseInt(pointsField.getText());
                String[] lines = pointsInput.getText().split("\n");
                double[] xPoints = new double[numPoints];
                double[] yPoints = new double[numPoints];
                for (int i = 0; i < numPoints; i++) {
                    String[] parts = lines[i].trim().split("\\s+");
                    if (parts.length != 2) throw new IllegalArgumentException("Invalid point format at line " + (i + 1));
                    xPoints[i] = Double.parseDouble(parts[0]);
                    yPoints[i] = Double.parseDouble(parts[1]);
                }
                double point = Double.parseDouble(xValueField.getText());
                String method = methodChoice.getValue();

                TextAreaOutputStream outStream = new TextAreaOutputStream(outputArea);
                System.setOut(new java.io.PrintStream(outStream));

                if (pointsChoice.getValue().equals("2 Points")) {
                    if (method.equalsIgnoreCase("forward")) {
                        differentiation.two_points_forward(numPoints, xPoints, yPoints, point);
                    } else {
                        differentiation.two_points_backward(numPoints, xPoints, yPoints, point);
                    }
                } else {
                    if (method.equalsIgnoreCase("forward")) {
                        differentiation.three_points_forward(numPoints, xPoints, yPoints, point);
                    } else if (method.equalsIgnoreCase("backward")) {
                        differentiation.three_points_backward(numPoints, xPoints, yPoints, point);
                    } else {
                        differentiation.central(numPoints, xPoints, yPoints, point);
                    }
                }
            } catch (Exception ex) {
                LOGGER.warning("Calculation error: " + ex.getMessage());
                showErrorDialog("Calculation Error", getUserFriendlyErrorMessage(ex));
                outputArea.setText("Error: " + getUserFriendlyErrorMessage(ex));
            } finally {
                progressIndicator.setVisible(false);
            }
        });

        addLabeledField(inputBox, "Output:", outputArea, "Differentiation results");
        inputBox.getChildren().add(calculateBtn);
        panel.getChildren().add(inputBox);
    }

    private VBox createInputBox() {
        VBox inputBox = new VBox(10);
        inputBox.setPadding(new Insets(15));
        inputBox.getStyleClass().add("input-box");
        return inputBox;
    }

    private TextField createTextField(String prompt) {
        TextField field = new TextField();
        field.setPromptText(prompt);
        field.getStyleClass().add("input-field");
        field.setTooltip(new Tooltip(prompt));
        addInputValidation(field);
        return field;
    }

    private TextArea createTextArea(String prompt) {
        TextArea area = new TextArea();
        area.setPromptText(prompt);
        area.getStyleClass().add("input-field");
        area.setPrefHeight(120);
        area.setTooltip(new Tooltip(prompt));
        return area;
    }

    private TextArea createOutputArea() {
        TextArea area = new TextArea();
        area.setEditable(false);
        area.getStyleClass().add("output-area");
        area.setPrefHeight(250);
        area.setTooltip(new Tooltip("Calculation results"));
        return area;
    }

    private ChoiceBox<String> createChoiceBox(String[] items, String tooltip) {
        ChoiceBox<String> choiceBox = new ChoiceBox<>();
        choiceBox.getItems().addAll(items);
        choiceBox.getStyleClass().add("input-field");
        choiceBox.setTooltip(new Tooltip(tooltip));
        return choiceBox;
    }

    private void addLabeledField(VBox box, String labelText, Region field, String tooltip) {
        Label label = new Label(labelText);
        label.getStyleClass().add("input-label");
        label.setTooltip(new Tooltip(tooltip));
        box.getChildren().addAll(label, field);
    }

    private void addInputValidation(TextField field) {
        field.textProperty().addListener((obs, oldValue, newValue) -> {
            try {
                if (!newValue.isEmpty() && field != null && !field.getPromptText().contains("function")) {
                    Double.parseDouble(newValue);
                    field.setStyle("-fx-border-color: #4CAF50;");
                }
            } catch (NumberFormatException e) {
                field.setStyle("-fx-border-color: #F44336;");
            }
        });
    }

    private void validateNonLinearInputs(TextField functionField, TextField initialField, TextField secondField, TextField toleranceField, TextField iterationsField, String method) {
        if (functionField.getText().isEmpty()) {
            throw new IllegalArgumentException("Function field cannot be empty");
        }
        if (initialField.getText().isEmpty()) {
            throw new IllegalArgumentException("Initial guess is required");
        }
        try {
            Double.parseDouble(initialField.getText());
        } catch (NumberFormatException e) {
            throw new IllegalArgumentException("Invalid initial guess");
        }
        if (method.equalsIgnoreCase("bisection") || method.equalsIgnoreCase("false position") || method.equalsIgnoreCase("secant")) {
            if (secondField.getText().isEmpty()) {
                throw new IllegalArgumentException("Second point is required");
            }
            try {
                Double.parseDouble(secondField.getText());
            } catch (NumberFormatException e) {
                throw new IllegalArgumentException("Invalid second point");
            }
        }
        if (method.equalsIgnoreCase("newton") || method.equalsIgnoreCase("halley") || method.equalsIgnoreCase("aitkens")) {
            if (toleranceField.getText().isEmpty()) {
                throw new IllegalArgumentException("Tolerance is required");
            }
            try {
                double tolerance = Double.parseDouble(toleranceField.getText());
                if (tolerance <= 0) throw new IllegalArgumentException("Tolerance must be positive");
            } catch (NumberFormatException e) {
                throw new IllegalArgumentException("Invalid tolerance value");
            }
        }
        if (iterationsField.getText().isEmpty()) {
            throw new IllegalArgumentException("Max iterations is required");
        }
        try {
            int iterations = Integer.parseInt(iterationsField.getText());
            if (iterations <= 0) throw new IllegalArgumentException("Iterations must be positive");
        } catch (NumberFormatException e) {
            throw new IllegalArgumentException("Invalid number of iterations");
        }
    }

    private void validateSystemInputs(ChoiceBox<String> typeChoice, TextField equationsField, TextArea coefficientsInput,
                                    TextField initialGuessesField, TextField toleranceField, TextField iterationsField, ChoiceBox<String> methodChoice) {
        if (typeChoice.getValue() == null) {
            throw new IllegalArgumentException("System type must be selected");
        }
        if (equationsField.getText().isEmpty()) {
            throw new IllegalArgumentException("Number of equations is required");
        }
        int n;
        try {
            n = Integer.parseInt(equationsField.getText());
            if (n <= 0) throw new IllegalArgumentException("Number of equations must be positive");
        } catch (NumberFormatException e) {
            throw new IllegalArgumentException("Invalid number of equations");
        }
        if (typeChoice.getValue().equals("Linear")) {
            String[] lines = coefficientsInput.getText().split("\n");
            if (lines.length < n + 1) throw new IllegalArgumentException("Insufficient input: provide " + n + " rows for A and 1 row for b");
            for (int i = 0; i < n; i++) {
                String[] parts = lines[i].trim().split("\\s+");
                if (parts.length != n) throw new IllegalArgumentException("Invalid row format for A at line " + (i + 1));
                for (String part : parts) {
                    try {
                        Double.parseDouble(part);
                    } catch (NumberFormatException e) {
                        throw new IllegalArgumentException("Invalid numerical value in A at line " + (i + 1));
                    }
                }
            }
            String[] bParts = lines[n].trim().split("\\s+");
            if (bParts.length != n) throw new IllegalArgumentException("Invalid format for b");
            for (String part : bParts) {
                try {
                    Double.parseDouble(part);
                } catch (NumberFormatException e) {
                    throw new IllegalArgumentException("Invalid numerical value in b");
                }
            }
        } else {
            String[] functions = coefficientsInput.getText().split("\n");
            if (functions.length != n) throw new IllegalArgumentException("Number of functions must match number of equations");
        }
        if (initialGuessesField.getText().isEmpty()) {
            throw new IllegalArgumentException("Initial guesses are required");
        }
        if (toleranceField.getText().isEmpty()) {
            throw new IllegalArgumentException("Tolerance is required");
        }
        try {
            double tolerance = Double.parseDouble(toleranceField.getText());
            if (tolerance <= 0) throw new IllegalArgumentException("Tolerance must be positive");
        } catch (NumberFormatException e) {
            throw new IllegalArgumentException("Invalid tolerance value");
        }
        if (iterationsField.getText().isEmpty()) {
            throw new IllegalArgumentException("Number of iterations is required");
        }
        try {
            int iterations = Integer.parseInt(iterationsField.getText());
            if (iterations <= 0) throw new IllegalArgumentException("Iterations must be positive");
        } catch (NumberFormatException e) {
            throw new IllegalArgumentException("Invalid number of iterations");
        }
        if (methodChoice.getValue() == null) {
            throw new IllegalArgumentException("Method must be selected");
        }
    }

    private void validateIntegrationInputs(TextField functionField, TextField upperField, TextField lowerField, TextField intervalsField, String method) {
        if (functionField.getText().isEmpty()) {
            throw new IllegalArgumentException("Function field cannot be empty");
        }
        if (upperField.getText().isEmpty() || lowerField.getText().isEmpty()) {
            throw new IllegalArgumentException("Upper and lower limits are required");
        }
        try {
            Double.parseDouble(upperField.getText());
            Double.parseDouble(lowerField.getText());
        } catch (NumberFormatException e) {
            throw new IllegalArgumentException("Invalid upper or lower limit");
        }
        if (method.equalsIgnoreCase("composite_trapezoidal") || method.equalsIgnoreCase("composite_simpson") || method.equalsIgnoreCase("romberg")) {
            if (intervalsField.getText().isEmpty()) {
                throw new IllegalArgumentException("Number of subintervals is required");
            }
            try {
                int intervals = Integer.parseInt(intervalsField.getText());
                if (intervals <= 0) throw new IllegalArgumentException("Subintervals must be positive");
            } catch (NumberFormatException e) {
                throw new IllegalArgumentException("Invalid number of subintervals");
            }
        }
    }

    private void validateCurveFittingInputs(TextField pointsField, TextArea pointsInput, TextField xValueField) {
        if (pointsField.getText().isEmpty()) {
            throw new IllegalArgumentException("Number of points is required");
        }
        try {
            int numPoints = Integer.parseInt(pointsField.getText());
            if (numPoints <= 0) throw new IllegalArgumentException("Number of points must be positive");
            String[] lines = pointsInput.getText().split("\n");
            if (lines.length < numPoints) throw new IllegalArgumentException("Insufficient points provided");
        } catch (NumberFormatException e) {
            throw new IllegalArgumentException("Invalid number of points");
        }
        if (xValueField.getText().isEmpty()) {
            throw new IllegalArgumentException("X value is required");
        }
        try {
            Double.parseDouble(xValueField.getText());
        } catch (NumberFormatException e) {
            throw new IllegalArgumentException("Invalid x value");
        }
    }

    private void validateLagrangeInputs(TextField pointsField, TextArea pointsInput, TextField xValueField) {
        if (pointsField.getText().isEmpty()) {
            throw new IllegalArgumentException("Number of points is required");
        }
        try {
            int n = Integer.parseInt(pointsField.getText());
            if (n <= 0) throw new IllegalArgumentException("Number of points must be positive");
            String[] lines = pointsInput.getText().split("\n");
            if (lines.length < n) throw new IllegalArgumentException("Insufficient points provided");
        } catch (NumberFormatException e) {
            throw new IllegalArgumentException("Invalid number of points");
        }
        if (xValueField.getText().isEmpty()) {
            throw new IllegalArgumentException("X value is required");
        }
        try {
            Double.parseDouble(xValueField.getText());
        } catch (NumberFormatException e) {
            throw new IllegalArgumentException("Invalid x value");
        }
    }

    private void validateDifferentiationInputs(ChoiceBox<String> pointsChoice, ChoiceBox<String> methodChoice,
                                            TextField pointsField, TextArea pointsInput, TextField xValueField) {
        if (pointsChoice.getValue() == null) {
            throw new IllegalArgumentException("Number of points must be selected");
        }
        if (methodChoice.getValue() == null) {
            throw new IllegalArgumentException("Method must be selected");
        }
        if (pointsField.getText().isEmpty()) {
            throw new IllegalArgumentException("Number of points is required");
        }
        try {
            int numPoints = Integer.parseInt(pointsField.getText());
            if (numPoints <= 0) throw new IllegalArgumentException("Number of points must be positive");
            String[] lines = pointsInput.getText().split("\n");
            if (lines.length < numPoints) throw new IllegalArgumentException("Insufficient points provided");
        } catch (NumberFormatException e) {
            throw new IllegalArgumentException("Invalid number of points");
        }
        if (xValueField.getText().isEmpty()) {
            throw new IllegalArgumentException("Point to differentiate is required");
        }
        try {
            Double.parseDouble(xValueField.getText());
        } catch (NumberFormatException e) {
            throw new IllegalArgumentException("Invalid point to differentiate");
        }
    }

    private String getUserFriendlyErrorMessage(Exception ex) {
        if (ex instanceof NumberFormatException) {
            return "Please enter valid numerical values.";
        } else if (ex instanceof IllegalArgumentException) {
            return ex.getMessage();
        } else {
            return "An unexpected error occurred: " + ex.getMessage() + "\nPlease check your inputs and try again.";
        }
    }

    private static class TextAreaOutputStream extends java.io.OutputStream {
        private TextArea textArea;

        public TextAreaOutputStream(TextArea textArea) {
            this.textArea = textArea;
        }

        @Override
        public void write(int b) {
            textArea.appendText(String.valueOf((char) b));
        }
    }
}