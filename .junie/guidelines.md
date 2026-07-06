# AI Developer Guidelines: Force Field X Development in IntelliJ IDEA with Junie

Welcome to the project! This document outlines the engineering standards, architecture, and IDE integrations required when modifying or expanding this numerical codebase.

## 1. Environment & IDE Standards
* **Target IDE**: IntelliJ IDEA (2026.1+). Always respect `.idea/` configuration files when present.
* **Java Version**: JDK 25 (LTS). Leverage modern features like Pattern Matching, Records for data structures, and Switch Expressions.
* **Build System**: Maven. External dependencies are managed through `pom.xml`.
* **External Libraries**:
    * See the parent 'pom.xml' and the those for each module.

## 2. Code Style & Formatting
* **Formatter**: Follow the integrated Google Java Style Guide. Run `Reformat Code` (`Ctrl+Alt+L` / `Cmd+Option+L`) before committing.
* **Naming Conventions**:
    * Math symbols must map clearly to code variable names (e.g., matrix $A$ becomes `matrixA`, vector $b$ becomes `vectorB`).
    * Use descriptive prefixes for accuracy variables (e.g., `tolerance` instead of `t`, `machineEpsilon` instead of `eps`).

## 3. Numerical Precision & Performance Guidelines
* **Floating-Point Types**: Use `double` for all standard calculations. Use `float` only if strictly required for memory/GPU constraints.
* **Object Allocation**: Minimize object creation inside tight loops (e.g., inside Runge-Kutta solvers or matrix iterations). Reuse mutable vector/matrix buffers where possible to prevent GC overhead.

## 4. Error Handling & Validation
* **Input Boundaries**: Validate numerical inputs at the public API boundary. Check for `Double.isNaN()`, `Double.isInfinite()`, and matrix dimension compatibility.
* **Exceptions**: Throw standard Java runtime exceptions:
    * `IllegalArgumentException` for mismatched dimensions or negative variances.
    * `ArithmeticException` for division by zero or non-convergent algorithms.

## 5. Testing Requirements
* **Framework**: JUnit 4.13.
* **Floating-Point Assertions**: Never use `assertEquals(expected, actual)`. Always use a tolerance threshold:
  ```java
  double tolerance = 1e-9;
  assertEquals(expected, actual, tolerance);
  ```

