
public class Main {

    public static void main(String[] args) {
        quickestDescentMethod();
    }

    //Метод Скорейшего Спуска
    static void quickestDescentMethod(){
        //Нужно найти s: посмотреть в методе золотого сечения замечание
        //|x2|, |x1| не стоит брать более 3 (экспонента будет возводиться в сумму их квадратов и выйдет за границы типа)
        Function func = new Function(5,1,-3,10e-6);
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

    //Первые производные: (+методы перегружены)
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

    //Вторые производные: (решил на листке, есть в телефоне фото) без перегрузки
    public double dx1dx1(double[] x){
        return 0;
    }
    public double dx2dx2(double[] x){
        return 0;
    }
    public double dx3dx3(double[] x){
        return 0;
    }
    public double dx1dx2(double[] x){
        return 0;
    }//=dx2dx1
    public double dx2dx3(double[] x){
        return 0;
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


}