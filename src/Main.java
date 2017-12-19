public class Main {

    public static void main(String...args) throws Exception {
        int nodeNum = 10;
        double lambda = 0.5;
        double intervalBottom = 0;
        double intervalUpper = 1;
        double step = (intervalUpper - intervalBottom) / nodeNum;
        double accuracy = Math.pow(step, 5);
        double[] nodes = new double [nodeNum + 1];

        for (int i = 0; i < nodeNum + 1; i++)
        {
            nodes[i] = intervalBottom + i * step;
        }

        System.out.println("Узлы:");
        printVector(nodes);

        System.out.println("Метод механических квадратур для интегрального уравнения Фредгольма:");
        System.out.println("Значения искомой функции в заданных узлах:");
        printVector(calcMechanicalQuadratureFredholm(nodes, step, lambda));
        System.out.println();

        System.out.println("Метод механических квадратур для интегрального уравнения Вольтерра:");
        System.out.println("Значения искомой функции в заданных узлах:");
        printVector(calcMechanicalQuadratureVolterra(nodes, step, lambda));
        System.out.println();

        System.out.println("Метод последовательных приближений для интегрального уравнения Фредгольма:");
        calcSequentialApproximateFredholm(nodes, step, lambda, accuracy);
        System.out.println();

        System.out.println("Метод последовательных приближений для интегрального уравнения Вольтерра");
        calcSequentialApproximateVolterra(nodes, step, lambda, accuracy);
        System.out.println();
    }

    private static void printVector(double[] vector) {
        for (int i = 0; i < vector.length; i++) {
            System.out.println(vector[i] + " ");
        }
    }

    private static double[] copyVector(double[] vector) {
        double[] newVect = new double[vector.length];
        for (int i = 0; i < vector.length; i++) {
            newVect[i] = vector[i];
        }
        return newVect;
    }

    private static double[] swapVectorElements(double[] vector, int from, int to) {
        double tmp = vector[from];
        vector[from] = vector[to];
        vector[to] = tmp;
        return vector;
    }

    private static double[] vectorDifference(double[] minuend, double[] subtrahend) {
        double[] difference = new double[minuend.length];
        for (int i = 0; i < difference.length; i++) {
            difference[i] = minuend[i] - subtrahend[i];
        }
        return difference;
    }

    private static double vectorNorm(double[] vector) {
        double max = vector[0];
        for (int i = 1; i < vector.length; i++) {
            if (vector[i] > max) {
                max = vector[i];
            }
        }
        return max;
    }

    private static double[] vectorPart(double[] vector, int indexFrom, int indexTo) throws Exception {
        double[] result;
        if (indexTo - indexFrom < 0 || indexFrom < 0 || indexTo < 0 || indexFrom >= vector.length || indexTo >= vector.length) {
            throw new Exception("Incorrect indices");
        }
        result = new double[indexTo - indexFrom + 1];
        for (int i = indexFrom; i <= indexTo; i++) {
            result[i - indexFrom] = vector[i];
        }
        return result;
    }

    private static double[][] copyMatrix(double[][] mtr) {
        double[][] newMtr = new double[mtr.length][mtr[0].length];
        for (int i = 0; i < mtr.length; i++) {
            for (int j = 0; j < mtr[i].length; j++) {
                newMtr[i][j] = mtr[i][j];
            }
        }
        return newMtr;
    }

    private static double[][] swapMatrixColumns(double[][] mtr, int from, int to) {
        double tmp;
        for (int i = 0; i < mtr.length; i++) {
            tmp = mtr[i][from];
            mtr[i][from] = mtr[i][to];
            mtr[i][to] = tmp;
        }
        return mtr;
    }

    private static double[] calcMechanicalQuadratureFredholm(double[] nodes, double step, double lambda) {
        double[] funcVect = new double[nodes.length];
        double[][] mtr = new double[nodes.length][nodes.length];
        double[] coefs = new double[nodes.length];

        for (int i = 0; i < coefs.length - 1; i++) {
            coefs[i] = step;
        }
        for (int i = 0; i < funcVect.length; i++) {
            funcVect[i] = calcFunction(nodes[i]);
        }
        for (int i = 0; i < mtr.length; i++) {
            for (int j = 0; j < mtr[i].length; j++) {
                mtr[i][j] = (i == j)
                        ? (1 - lambda * coefs[j] * calcKernel(nodes[i], nodes[j]))
                        : (-lambda * coefs[j] * calcKernel(nodes[i], nodes[j]));
            }
        }

        return gauss(mtr, funcVect);
    }

    private static double calcFunction(double x) {
        return 2 - Math.sin(Math.PI * x);
    }

    private static double calcKernel(double x, double s) {
        return 1 / (2 + Math.sin(Math.PI * (x + s)));
    }

    private static double[] gauss(double[][] mtrA, double[] vectB) {
        double[][] a = copyMatrix(mtrA);
        double[] b = copyVector(vectB);
        double[] x = new double[b.length];
        double[] xIndexes = new double[b.length];
        double max;
        int maxK;
        int index;

        for (int i = 0; i < x.length; i++) {
            xIndexes[i] = (new Integer(i)).doubleValue();
        }

        for (int k = 0; k < x.length; k++) {
            max = a[k][k];
            maxK = k;
            for (int i = k + 1; i < x.length; i++) {
                if (Math.abs(max) < Math.abs(a[k][i])) {
                    max = a[k][i];
                    maxK = i;
                }
            }
            if (maxK != k) {
                a = swapMatrixColumns(a, k, maxK);
                xIndexes = swapVectorElements(xIndexes, k, maxK);
            }
            for (int j = k; j < x.length; j++) {
                a[k][j] /= max;
            }
            b[k] /= max;
            for (int i = k + 1; i < x.length; i++) {
                for (int j = k + 1; j < x.length; j++) {
                    a[i][j] -= a[i][k] * a[k][j];
                }
                b[i] -= a[i][k] * b[k];
                a[i][k] = 0.0;
            }
        }

        for (int i = x.length - 1; i >= 0; i--) {
            index = (new Double(xIndexes[i])).intValue();
            x[index] = b[i];
            for (int j = i + 1; j < x.length; j++) {
                x[index] -= a[i][j] * x[(new Double(xIndexes[j])).intValue()];
            }
        }

        return x;
    }

    private static double[] calcMechanicalQuadratureVolterra(double[] nodes, double step, double lambda) {
        double[] result = new double[nodes.length];
        double sum;

        for (int i = 0; i < nodes.length; i++) {
            sum = 0;
            for (int j = 0; j < i; j++) {
                sum += calcKernel(nodes[i], nodes[j]) * step * result[j];
            }
            result[i] = (calcFunction(nodes[i]) + lambda * sum) / (1 - lambda * step * calcKernel(nodes[i], nodes[i]));
        }
        return result;
    }

    private static void calcSequentialApproximateFredholm(double[] nodes, double step, double lambda, double accuracy) {
        double[] currentResult;
        double[] nextResult = new double[nodes.length];
        int counter = 0;

        for (int i = 0; i < nextResult.length; i++) {
            nextResult[i] = calcFunction(nodes[i]);
        }

        do {
            currentResult = copyVector(nextResult);
            for (int i = 0; i < nextResult.length; i++) {
                nextResult[i] = lambda * calcIntegralLeftRectangles(i, nodes, currentResult, step) + calcFunction(nodes[i]);
            }
            counter++;
        } while (vectorNorm(vectorDifference(nextResult, currentResult)) >= accuracy);
        System.out.println("Количество итераций: " + counter);
        System.out.println("Значения искомой функции в заданных узлах:");
        printVector(nextResult);
    }

    private static double calcIntegralLeftRectangles(int index, double[] nodes, double[] solution, double step) {
        double sum = 0;
        for (int i = 0; i < nodes.length - 1; i++) {
            sum += calcKernel(nodes[index], nodes[i]) * solution[i];
        }
        return step * sum;
    }

    private static void calcSequentialApproximateVolterra(double[] nodes, double step, double lambda, double accuracy) throws Exception {
        double[] currentResult;
        double[] nextResult = new double[nodes.length];
        int counter = 0;

        for (int i = 0; i < nextResult.length; i++) {
            nextResult[i] = calcFunction(nodes[i]);
        }
        do {
            currentResult = copyVector(nextResult);
            for (int i = 0; i < nextResult.length; i++) {
                nextResult[i] = lambda * calcIntegralMiddleRectanglesVolterra(i, nodes, currentResult, step) + calcFunction(nodes[i]);
            }
            counter++;
        } while (vectorNorm(vectorDifference(nextResult, currentResult)) >= accuracy);
        System.out.println("Количество итераций: " + counter);
        System.out.println("Значения искомой функции в заданных узлах:");
        printVector(nextResult);
    }

    private static double calcIntegralMiddleRectanglesVolterra(int index, double[] nodes, double[] solution, double step) throws Exception {
        double[] solutionMiddle = new double[solution.length];
        for (int i = 0; i < solution.length; i++) {
            solutionMiddle[i] = interpolateLagrange(nodes, solution, nodes[i] + step / 2);
        }
        return calcIntegralMiddleRectangles(index, vectorPart(nodes, 0, index), vectorPart(solutionMiddle, 0, index), step);
    }

    private static double interpolateLagrange(double[] nodes, double[] values, double x) {
        double result = 0;
        double tmp;
        for (int i = 0; i < nodes.length; i++) {
            tmp = values[i];
            for (int j = 0; j < nodes.length; j++) {
                tmp *= (i != j) ? ((x - nodes[j]) / (nodes[i] - nodes[j])) : 1;
            }
            result += tmp;
        }
        return result;
    }

    private static double calcIntegralMiddleRectangles(int index, double[] nodes, double[] solution, double step) {
        double sum = 0;
        for (int i = 0; i < nodes.length - 1; i++) {
            sum += calcKernel(nodes[index], nodes[i] + step / 2) * interpolateLagrange(nodes, solution, nodes[i] + step / 2);
        }
        return step * sum;
    }
}
