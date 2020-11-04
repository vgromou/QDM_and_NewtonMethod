import java.util.Arrays;

public class Main {

    public static void main(String[] args) {
        quickestDescentMethod();
        System.out.println("\n\n\n");
        newtonMethod();
    }

    //Метод Скорейшего Спуска
    static void quickestDescentMethod(){
        //Нужно найти s: посмотреть в методе золотого сечения замечание
        //|x2|, |x1| не стоит брать более 3 (экспонента будет возводиться в сумму их квадратов и выйдет за границы типа)
        //Деления на ноль нет
        Function func = new Function(1,2,1,10e-6);
        double sk;
        System.out.println("\t\t\t\t\t\tМЕТОД СКОРЕЙШЕГО СПУСКА\nДанные по каждой итерации:");
        System.out.println("#####   x                                     f(x)        |f'(x)|");
        System.out.printf("%5d   [% 3.6f, % 3.6f, % 3.6f]    % 3.6f    %3.6f\n", 0, func.x[0], func.x[1], func.x[2],
                func.f(func.x[0],func.x[1],func.x[2]), func.dxNorm());
        for(int i = 1; func.gLength > func.eps; i++){
            sk = func.findArgMin(0.1); //h = шаг
            func.setX(sk);
            func.g[0] = func.dx1(func.x);
            func.g[1] = func.dx2(func.x);
            func.g[2] = func.dx3(func.x);

            func.gLength = func.nG();
            System.out.printf("%5d   [% 3.6f, % 3.6f, % 3.6f]    % 3.6f    %3.6f\n", i, func.x[0], func.x[1], func.x[2],
                    func.f(func.x[0],func.x[1],func.x[2]), func.gLength);
        }
        System.out.println("______________________________________________________________________________");
        System.out.println("Результат:");
        System.out.println("   xMin: [" + func.x[0] + ", " + func.x[1] + ", " + func.x[2] + "]");
        System.out.println("f(xMin): " + func.f(func.x));

    }

    //Метод Ньютона
    static void newtonMethod(){
        //Не любой вектор пройдет. Иногда, при поиске обратной матрицы, определитель будет равен нулю
        //т.е. произойдет деление на ноль и получится бесконечность. На такой случай введена проверка равенства
        //определителя нулю (с точностью заданной эпсилон), и если он нулю равен, вылетит арифметическая ошибка
        Function func = new Function(1,2,1,10e-10);
        System.out.println("\t\t\t\t\t\t\t\tМЕТОД НЬЮТОНА\nДанные по каждой итерации:");
        System.out.println("#####   x                                     f(x)        |f'(x)|");
        System.out.printf("%5d   [% 3.6f, % 3.6f, % 3.6f]    % 3.6f    %3.6f\n", 0, func.x[0], func.x[1], func.x[2],
                func.f(func.x[0],func.x[1],func.x[2]), func.dxNorm());
        double[] dx = func.getFdx();
        double[][] dxdx = Matrix.getMatrix(func, func.x);
        for(int i = 1; func.dxNorm() > func.eps; i++){
            double[][] inversed = Matrix.inverse(dxdx);
            func.setX(Matrix.multiplyByVector(inversed, dx));
            dx = func.getFdx();
            dxdx = Matrix.getMatrix(func, func.x);
            System.out.printf("%5d   [% 3.6f, % 3.6f, % 3.6f]    % 3.6f    %3.6f\n", i, func.x[0], func.x[1], func.x[2],
                    func.f(func.x[0],func.x[1],func.x[2]), func.dxNorm());
        }
        System.out.println("______________________________________________________________________________");
        System.out.println("Результат:");
        System.out.println("   xMin: [" + func.x[0] + ", " + func.x[1] + ", " + func.x[2] + "]");
        System.out.println("f(xMin): " + func.f(func.x));
    }
}

class Function{
    //function: e^(x1^2+x2^2)+ln(4+x2^2+2*x3^2)
    //EPS = 10e-6 для метода скорейшего спуска
    //EPS = 10e-10 для метода Ньютона
    double eps;
    public double[] x;
    public double[] g = new double[3];
    public double gLength;

    public Function(double x1, double x2, double x3, double eps){
        x = new double[]{x1,x2,x3};
        g[0] = dx1(x);
        g[1] = dx2(x);
        g[2] = dx3(x);
        gLength = nG(g[0],g[1],g[2]);
        this.eps = eps;
    }

    //Функция: (+перегружен)
    public double f(double x1, double x2, double x3) {
        return (Math.exp(x1 * x1 + x2 * x2) + Math.log(4 + x2 * x2 + 2 * x3 * x3));
    }
    public double f(double[] x){
        return f(x[0], x[1], x[2]);
    }

    //Первые производные: (+перегружены)
    public double dx1(double x1, double x2) {
        return (2 * x1 * Math.exp(x1 * x1 + x2 * x2));
    }
    public double dx2(double x1, double x2, double x3){
        return (2 * x2 * Math.exp(x1 * x1 + x2 * x2) + 2 * x2 / (4 + x2 * x2 + 2 * x3 * x3));
    }
    public double dx3(double x2, double x3){
        return (4 * x3 / (4 + x2 * x2 + 2 * x3 * x3));
    }
    public double dx1(double[] x){
        return dx1(x[0],x[1]);
    }
    public double dx2(double[] x){
        return dx2(x[0],x[1],x[2]);
    }
    public double dx3(double[] x){
        return dx3(x[1],x[2]);
    }

    //Вторые производные:
    //x[0] = x1; x[1] = x2; x[2] = x3
    public double dx1dx1(double[] x){
        return (2 * Math.exp(x[0] * x[0] + x[1] * x[1]) + 4 * x[0] * x[0] * Math.exp(x[0] * x[0] + x[1] * x[1]));
    }
    public double dx2dx2(double[] x){
        double e = Math.exp(x[0] * x[0] + x[1] * x[1]);
        double del = (4 + x[1] * x[1] + 2 * x[2] * x[2]);
        return (2 * e + 4 * x[1] * x[1] * e + (8 - 2 * x[1] * x[1] + 4 * x[2] * x[2])/(del * del));
    }
    public double dx3dx3(double[] x){
        double del = (4 + x[1] * x[1] + 2 * x[2] * x[2]);
        return ((16 + 4 * x[1] * x[1] - 8 * x[2] * x[2])/(del * del));
    }
    public double dx1dx2(double[] x){
        return (4 * x[0] * x[1] + Math.exp(x[0] * x[0] + x[1] * x[1]));
    }//=dx2dx1
    public double dx2dx3(double[] x){
        double del = (4 + x[1] * x[1] + 2 * x[2] * x[2]);
        return (-8 * x[1] * x[2] / (del * del));
    }//=dx3dx2
    public double dx3dx1(double[] x){
        return 0;
    }//=dx1dx3
    public double dx2dx1(double[] x){ return dx1dx2(x);}
    public double dx3dx2(double[] x){ return dx2dx3(x);}
    public double dx1dx3(double[] x){ return dx3dx1(x);}

    //Норма g: (+перегружен)
    public double nG(double gk1, double gk2, double gk3) {
        return (Math.sqrt(gk1 * gk1 + gk2 * gk2 + gk3 * gk3));
    }
    public double nG() {
        return nG(g[0], g[1], g[2]);
    }

    public void setX(double s){
        x[0] = x[0] - s * g[0];
        x[1] = x[1] - s * g[1];
        x[2] = x[2] - s * g[2];
    }//Для Метода Скорейшего Спуска
    public void setX(double[] v){
        x[0] = x[0] - v[0];
        x[1] = x[1] - v[1];
        x[2] = x[2] - v[2];
    }//Для Метода Ньютона
    public double phi(double s){
        double x1 = x[0] - s * g[0];
        double x2 = x[1] - s * g[1];
        double x3 = x[2] - s * g[2];
        return f(x1, x2, x3);
    }

    //Поиск Sk
    public double findArgMin(double h){
        double a = 0;
        double b = findB(h);

        double delta = ((b-a)*(3 - Math.sqrt(5)))/2;
        double x = a + delta;
        double y = b - delta;

        while(b - a >= 2 * eps){
            if(phi(x) > phi(y)){
                a = x;
                x = y;
                y = b + a - y;
            }
            else{
                b = y;
                y = x;
                x = a + b - y;
            }
        }
        return (a+b)/2;
    }
    //Поиск правой границы с шагом h
    double findB(double h){
        //h - шаг
        double phiA = phi(0);
        double phiH = phi(h);
        while(phiA >= phiH){
            h = 2 * h;
            phiH = phi(h);
        }
        return h;
    }

    // |f(x)|
    public double norm(){return Math.sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);}
    // |f'(x)|
    public double dxNorm(){
        double d1 = this.dx1(x);
        double d2 = this.dx2(x);
        double d3 = this.dx3(x);
        return Math.sqrt(d1 * d1 + d2 * d2 + d3 * d3);
    }
    public double[] getFdx(){
        return new double[]{this.dx1(x), this.dx2(x), this.dx3(x)};
    }
}

abstract class Matrix{

    public static double[][] inverse(double[][] matrix){
        double determinant = findDeterminant(matrix);
        double[][] complement = getComplementMatrix(matrix);
        double[][] transposed = transpose(complement);
        if(Math.abs(determinant) < 10e-10) throw new ArithmeticException();
        return multiplyByNumber(transposed,1/determinant);
    }

    public static double[][] transpose(double[][] matrix){
        double[][] res = new double[matrix.length][matrix[0].length];
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[0].length; j++) {
                res[i][j] = matrix[j][i];
            }
        }
        return res;
    }

    public static double[][] multiplyByNumber(double[][] matrix, double n){
        double[][] res = new double[matrix.length][matrix[0].length];
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[0].length; j++) {
                res[i][j] = matrix[i][j] * n;
            }
        }
        return res;
    }

    public static double[] multiplyByVector(double[][] matrix, double[] vector){
        double[] res = new double[vector.length];
        for (int i = 0; i < matrix.length; i++) {
            res[i] = matrix[i][0] * vector[0];
            for (int j = 1; j < matrix[0].length; j++) {
                res[i] += matrix[i][j] * vector[j];
            }
        }
        return res;
    }

    public static double[][] getMatrix(Function f, double[] x){
        return new double[][] {{f.dx1dx1(x), f.dx1dx2(x), f.dx1dx3(x)},
                               {f.dx2dx1(x), f.dx2dx2(x), f.dx2dx3(x)},
                               {f.dx3dx1(x), f.dx3dx2(x), f.dx3dx3(x)}};
    }

    public static double[][] getComplementMatrix(double[][] matrix){
        double[][] res = new double[matrix.length][matrix[0].length];
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[0].length; j++) {
                int sign;
                if((i + j) %2 == 0) sign = 1;
                else sign = -1;
                res[i][j] = sign * findDeterminant(getMinor(matrix,i,j));
            }
        }
        return res;
    }

    public static double[][] getMinor(double[][] matrix, int row, int column){
        double[][] res = new double[matrix.length - 1][matrix[0].length - 1];
        for (int i = 0, k = 0; i < matrix.length; i++) {
            if(i == row) continue;
            for (int j = 0, m = 0; j < matrix[0].length; j++) {
                if(j == column) continue;
                res[k][m] = matrix[i][j];
                m++;
            }
            k++;
        }
        return res;
    }

    public static double findDeterminant(double[][] matrix){
        double res = 0;
        if(matrix.length < 3) return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
        for (int i = 0; i < matrix.length; i++) {
            double[][] temp = getMinor(matrix, 0, i);
            double detMinor = findDeterminant(temp);
            if (i % 2 == 0) res += matrix[0][i] * detMinor;
            else res -= matrix[0][i] * detMinor;
        }
        return res;
    }
}