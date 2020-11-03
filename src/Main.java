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
        Function func = new Function(1,3,10,10e-6);
        double sk;
        System.out.println("\t\t\t\t\t\tМЕТОД СКОРЕЙШЕГО СПУСКА\nДанные по каждой итерации:");
        System.out.println("#####   x                                     f(x)        |f'(x)|");
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
        //Нужны какие-то определенные значения x, иначе та же проблема, что в методе скорейшего спуска
        //но уже в середине
        Function func = new Function(0,0.1,2,10e-10);
        double[] funcDx = func.g; //f'(x)
        Matrix matrix = new Matrix();
        matrix.calcMatrix(func, func.x); //f''(x) (считает значение второй производной в точке x)
        System.out.println("\t\t\t\t\t\t\tМЕТОД НЬЮТОНА\nДанные по каждой итерации:");
        System.out.println("#####   x                                     f(x)        |f'(x)|");
        for(int i = 1; func.norm() > func.eps; i++){
            double[][] inverseMatrix = matrix.inverseMatrix();
            double[] v = Matrix.multiplyByVector(inverseMatrix, funcDx);
            func.setX(v);
            funcDx = new double[]{func.dx1(func.x), func.dx2(func.x), func.dx3(func.x)};
            matrix.calcMatrix(func, func.x);
            System.out.printf("%5d   [% 3.6f, % 3.6f, % 3.6f]    % 3.6f    %3.6f\n", i, func.x[0], func.x[1], func.x[2],
                    func.f(func.x[0],func.x[1],func.x[2]), func.norm());
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
    }
    public void setX(double[] v){
        x[0] = x[0] - v[0];
        x[1] = x[1] - v[1];
        x[2] = x[2] - v[2];
    }
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

    // |f'(x)|
    public double norm(){return nG(x[0],x[1],x[2]);}
}

class Matrix{
    //Класс подходит только для этого задания (3*3, вторые частные производные Function)
    double[][] matrix;
    double detMatrix;
    double[][] inverse;
    double[][] trans;

    public Matrix(){
        matrix = new double[3][3];
        inverse = new double[3][3];
    }

    public void calcMatrix(Function f, double[] x){
        this.matrix = new double[][] {{f.dx1dx1(x), f.dx1dx2(x), f.dx3dx1(x)},
                                      {f.dx1dx2(x), f.dx2dx2(x), f.dx2dx3(x)},
                                      {f.dx3dx1(x), f.dx2dx3(x), f.dx3dx3(x)}};
    }

    void calcDeterminant(){
        double[][] m = this.matrix;
        this.detMatrix = m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1])-
                m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0]) + m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]);
    }

    public double[][] inverseMatrix(){
        this.transpose();
        this.calcDeterminant();
        this.inverse = multiplyByNumber(this.trans, 1/this.detMatrix);
        return this.inverse;
    }

    public static double[][] multiplyByNumber(double[][] m, double k){
        m[0][0] = m[0][0] * k;
        m[0][1] = m[0][1] * k;
        m[0][2] = m[0][2] * k;
        m[1][0] = m[1][0] * k;
        m[1][1] = m[1][1] * k;
        m[1][2] = m[1][2] * k;
        m[2][0] = m[2][0] * k;
        m[2][1] = m[2][1] * k;
        m[2][2] = m[2][2] * k;
        return m;
    }

    void transpose(){
        double[][] m = this.matrix;
        this.trans = new double[3][3];
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                this.trans[i][j] = m[j][i];
            }
        }
    }

    public static double[] multiplyByVector(double[][] m, double[] v){
        double row1 = m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2];
        double row2 = m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2];
        double row3 = m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2];
        return new double[]{row1, row2, row3};
    }

    @Override
    public String toString() {
        return "Matrix{" +
                "matrix=" + Arrays.deepToString(matrix) +
                '}';
    }
}